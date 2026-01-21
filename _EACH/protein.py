import os
from Bio.Seq import Seq
from Bio import SeqUtils
from modules.signal import getSignalProteome
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint
# Config variables
N_DIGITS_GLOBAL = 3

class Protein:
    
    # ---------------------------- Static Global Variables  ------------------------
    # ------------------------------------------------------------------------------
    AMINO_ACID_WEIGHTS = {
        #Canonical 20 amino acids
        'A': 89.1, 'R': 174.2, 'N': 132.1, 'D': 133.1, 'C': 121.2, 'E': 147.1,
        'Q': 146.2, 'G': 75.1, 'H': 155.2, 'I': 131.2, 'L': 131.2, 'K': 146.2,
        'M': 149.2, 'F': 165.2, 'P': 115.1, 'S': 105.1, 'T': 119.1, 'W': 204.2,
        'Y': 181.2, 'V': 117.1,
        # Non canonical
        'U': 168.1, 'O': 255.3,
    }
    # Ambiguous codes (average weights)
    AMINO_ACID_WEIGHTS.update({
        'B': (AMINO_ACID_WEIGHTS['D'] + AMINO_ACID_WEIGHTS['N']) / 2,  # Aspartic acid or Asparagine
        'Z': (AMINO_ACID_WEIGHTS['E'] + AMINO_ACID_WEIGHTS['Q']) / 2,  # Glutamic acid or Glutamine
        'J': (AMINO_ACID_WEIGHTS['I'] + AMINO_ACID_WEIGHTS['L']) / 2,  # Isoleucine or Leucine
        'X': sum(AMINO_ACID_WEIGHTS.values()) / len(AMINO_ACID_WEIGHTS)  # Any amino acid
    })
    #  	Kyte-Doolittle hydrophobicity scale

    AMINO_ACID_HYDROPHOBICITY = {
        "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5, "M": 1.9, 
        "A": 1.8, "G": -0.4, "T": -0.7, "S": -0.8, "W": -0.9, 
        "Y": -1.3, "P": -1.6, "H": -3.2, "E": -3.5, "Q": -3.5, "D": -3.5,
        "N": -3.5, "K": -3.9, "R": -4.5
    }
    
    # Ambiguous codes (average hydrophobicity)
    AMINO_ACID_HYDROPHOBICITY.update({
        'B': (AMINO_ACID_HYDROPHOBICITY['D'] + AMINO_ACID_HYDROPHOBICITY['N']) / 2,  # Aspartic acid or Asparagine
        'Z': (AMINO_ACID_HYDROPHOBICITY['E'] + AMINO_ACID_HYDROPHOBICITY['Q']) / 2,  # Glutamic acid or Glutamine
        'J': (AMINO_ACID_HYDROPHOBICITY['I'] + AMINO_ACID_HYDROPHOBICITY['L']) / 2,  # Isoleucine or Leucine
        'X': sum(AMINO_ACID_HYDROPHOBICITY.values()) / len(AMINO_ACID_HYDROPHOBICITY)  # Any amino acid
    })

    # Amino acids polarity weights
    AMINO_ACID_POLARITY = {
        "R": 0, "K": 0, "D": 0.9, "E": 0.9, "N": 1, "Q": 1, "S": 1, "T": 1,
        "A": 0, "G": 0, "P": 0, "C": 0.55, "V": 0, "I": 0, "L": 0,
        "F": 0, "Y": 0.65, "W": 0.65, "H": 0.7, "M": 0
    }
    
    childClasses = {}
    
    proteomeIdLookupTable = {
        9606: "UP000005640",  # Homo sapiens
        3702: "UP000006548",  # Arabidopsis thaliana
    }
    
    masterProteomeID = None
    
    # ------------------------------------------------------------------------------
    # ---------------------------- Static Global Variables  ------------------------


    # ------------------------------- Static methods ------------------------------
    # ------------------------------------------------------------------------------
    @staticmethod
    def getAllProteins():
        """
        Return all instantiated Protein objects (registry store).

        :return: List of Protein instances currently tracked in `childClasses`.

        """
        return list(Protein.childClasses.values())
    
    @staticmethod
    def deleteAllProteins():
        """
        Clear the global registry of Protein instances.

        :return: None. Empties `childClasses`.

        """
        Protein.childClasses = {}
        Protein.masterProteomeID = None
    
    
    @staticmethod
    def depleteProteinsByWeight(minWeight=None, maxWeight=None):
        """
        Module helper (molecular_weight_cutoff): zero abundance outside weight bounds.

        :param minWeight: Minimum weight in kDa to keep (None disables lower bound).
        :param maxWeight: Maximum weight in kDa to keep (None disables upper bound).
        :return: None. Mutates Protein instances and logs modifications.

        """
        for protein in Protein.childClasses.values():
            if minWeight is not None and protein.get_weight() < minWeight:
                protein.abundance = 0.0
                protein.modifications.append(f"Depleted due to weight < {minWeight} kDa")
            if maxWeight is not None and protein.get_weight() > maxWeight:
                protein.abundance = 0.0
                protein.modifications.append(f"Depleted due to weight > {maxWeight} kDa")
    
                
    @staticmethod
    def fractionateProteinsByMolecularWeight(keepInsideOutsideSelection="inside", minWeight=None, maxWeight=None):
        """
        Module helper (SDS_page_fractionation): deplete proteins inside or outside a weight window.

        :param keepInsideOutsideSelection: "inside" to remove inside range, "outside" to remove outside.
        :param minWeight: Lower bound of weight window (None treated as -inf).
        :param maxWeight: Upper bound of weight window (None treated as +inf).
        :return: None. Mutates Protein instances and logs modifications.

        """
        if minWeight is None:
            minWeight = float('inf')
        if maxWeight is None:
            maxWeight = -float('inf')
        if keepInsideOutsideSelection == "outside":    
            for protein in Protein.childClasses.values():
                if minWeight <= protein.get_weight() <= maxWeight:
                    protein.abundance = 0.0
                    protein.modifications.append(f"Removed due to weight outside range {minWeight}-{maxWeight} kDa")
        elif keepInsideOutsideSelection == "inside":
            for protein in Protein.childClasses.values():
                if protein.get_weight() < minWeight:
                    protein.abundance = 0.0
                    protein.modifications.append(f"Removed due to weight < {minWeight} kDa")
                if protein.get_weight() > maxWeight:
                    protein.abundance = 0.0
                    protein.modifications.append(f"Removed due to weight > {maxWeight} kDa")
                
    @staticmethod
    def fractionateProteinsByIsoelectricPoint(keepInsideOutsideSelection="inside", minPI=None, maxPI=None):
        """
        Module helper (isoelectric_focussing): deplete proteins inside or outside a pI window.

        :param keepInsideOutsideSelection: "inside" to remove inside range, "outside" to remove outside.
        :param minPI: Lower bound of pI window (None treated as -inf).
        :param maxPI: Upper bound of pI window (None treated as +inf).
        :return: None. Mutates Protein instances and logs modifications.

        """
        if maxPI is None:
            maxPI = float('inf')
        if minPI is None:
            minPI = -float('inf')
        if keepInsideOutsideSelection == "outside":    
            for protein in Protein.childClasses.values():
                if minPI <= protein.isoelectric_point <= maxPI:
                    protein.abundance = 0.0
                    protein.modifications.append(f"Removed due to pI outside range {minPI}-{maxPI}")
        elif keepInsideOutsideSelection == "inside":
            for protein in Protein.childClasses.values():
                if protein.isoelectric_point < minPI:
                    protein.abundance = 0.0
                    protein.modifications.append(f"Removed due to pI < {minPI}")
                if protein.isoelectric_point > maxPI:
                    protein.abundance = 0.0
                    protein.modifications.append(f"Removed due to pI > {maxPI}")
                    
        
    @staticmethod
    def fractionateProteinsByHydrophobicity(keepInsideOutsideSelection="inside", minHydrophobicity=None, maxHydrophobicity=None):
        """
        Module helper (reversed_phase_chromatography): deplete proteins inside or outside a hydrophobicity window.

        :param keepInsideOutsideSelection: "inside" to remove inside range, "outside" to remove outside.
        :param minHydrophobicity: Lower bound of hydrophobicity window (None treated as -inf).
        :param maxHydrophobicity: Upper bound of hydrophobicity window (None treated as +inf).
        :return: None. Mutates Protein instances and logs modifications.

        """
        if maxHydrophobicity is None:
            maxHydrophobicity = float('inf')
        if minHydrophobicity is None:
            minHydrophobicity = -float('inf')
        
        if keepInsideOutsideSelection == "outside":    
            for protein in Protein.childClasses.values():
                if minHydrophobicity <= protein.hydrophobicity <= maxHydrophobicity:
                    protein.abundance = 0.0
                    protein.modifications.append(f"Removed due to hydrophobicity outside range {minHydrophobicity}-{maxHydrophobicity}")
        elif keepInsideOutsideSelection == "inside":
            for protein in Protein.childClasses.values():
                if protein.hydrophobicity < minHydrophobicity:
                    protein.abundance = 0.0
                    protein.modifications.append(f"Removed due to hydrophobicity < {minHydrophobicity}")
                if protein.hydrophobicity > maxHydrophobicity:
                    protein.abundance = 0.0
                    protein.modifications.append(f"Removed due to hydrophobicity > {maxHydrophobicity}")
                    
        
    @staticmethod
    def signalPeptideCleavage():
        """
        Module helper (signal_peptide_removal): cleave signal peptides using configured proteome map.

        :return: None. Mutates Protein sequences, updates attributes, logs modifications.
        """
        signalPeptides = getSignalProteome(proteomeId=Protein.masterProteomeID)
        for proteinName,startEnd in signalPeptides["protein_name"].items():
            proteinClass = Protein.childClasses.get(proteinName)
            if proteinClass:
                if proteinClass.get_abundance() > 0.01:
                    print(proteinClass, startEnd)
                proteinClass.__cleaveSequence(startEnd)
    
    @staticmethod
    def proteinImmunoaffinityDepletion(depletionDictionary):
        """
        Module helper (affinity_depletion): apply abundance depletion by entry name.

        :param depletionDictionary: Mapping entryName -> depletion fraction (0-1).
        :return: None. Mutates Protein instances and logs modifications.
        """
        for entryNameKey, proteinDepletionValue in depletionDictionary.items():
            if entryNameKey in Protein.childClasses:
                Protein.childClasses[entryNameKey].__depleteAbundance(proteinDepletionValue)
                
    @staticmethod  
    def saveProteinsAsFasta(filepath, proteinList=None):
        """
        Write proteins to a FASTA file (utility, not module-specific).

        :param filepath: Output path without extension.
        :param proteinList: Iterable of Protein or dict of entryName->Protein (defaults to all tracked proteins).
        :return: None. Writes `<filepath>.fasta`.
        """
        if proteinList is None:
            proteinList = Protein.childClasses
        if isinstance(proteinList, dict):
            proteinList = list(proteinList.values())

       # print(proteinList)
        with open(f"{filepath}.fasta", 'w') as fasta_file:
            for protein in proteinList:
                fasta_file.write(protein.get_fasta() + "\n")
                
    # ------------------------------------------------------------------------------
    # ------------------------------- Static methods  ------------------------------

    
    # ------------------------------ Internal methods ------------------------------
    # ------------------------------------------------------------------------------
    def __init__(self, header, sequence):
        """
        Initialize a Protein instance and register it in the global registry.

        Internal note: this constructor registers the instance in `childClasses`; avoid manual edits.

        :param header: Raw protein header (e.g., FASTA header).
        :param sequence: Amino-acid sequence (string) converted to Bio.Seq.Seq.
        :return: None.

        """
        
        self.header = header
        self.sequence = Seq(sequence)
        # Set Weight and other sequence-dependent attributes
        self.setSequenceDependentAttributes()
        self.__processHeader()
        self.modifications = []
        Protein.childClasses[self.entryName] = self
        
    def __processHeader(self):
        """
        INTERNAL: parse header tags and populate metadata fields.

        Warning: class-specific; callers should not mutate header parsing.

        Parse the instance header (self.header) and populate header-related attributes.
        Expected header format: "<db>|<accession>|<entryName> <TAG=VALUE> ...".
        This method:
        - Splits the header on '|' to set self.db and self.accession, then splits the remainder on spaces.
        - Sets self.entryName to the first token of the remainder.
        - Extracts tagged fields from tokens of the form KEY=VALUE and assigns:
            - OX=<int> -> self.organismID (int); if 9606, also sets self.organism = "Homo Sapiens"
            - GN=<str>  -> self.geneName (str)
            - PE=<int>  -> self.proteinExistence (int)
            - SV=<int>  -> self.sequenceVersion (int)
            - AB=<float>-> self.abundance (float; defaults to 0.0)
        - Leaves self.proteinName as None (not parsed here).
        - Ignores unrecognized tokens. Note: numeric conversions may raise ValueError on invalid input.
        
        :return: None. Sets db, accession, entryName, organism info, gene name, PE, SV, abundance.
        """
        self.db, self.accession, remainder = self.header.split('|',2)
        
        remainder = remainder.split(' ')

        self.entryName: str = None
        self.proteinName: str = None
        self.organism: str = None
        self.organismID: int = None
        self.geneName: str = None
        self.proteinExistence: int = None
        self.sequenceVersion: int = None
        self.abundance: float = 0.0
        self.entryName: str = remainder[0]    

        for item in remainder:
            if item.startswith('OX='):
                self.organismID = int(item[3:])
                if self.organismID == 9606:
                    self.organism = "Homo Sapiens"
            elif item.startswith('GN='):
                self.geneName = item[3:]
            elif item.startswith('PE='):
                self.proteinExistence = int(item[3:])
            elif item.startswith('SV='):
                self.sequenceVersion = int(item[3:])
            elif item.startswith('AB='):
                self.abundance = float(item[3:])
        
        if self.organismID and Protein.masterProteomeID is None:
            Protein.masterProteomeID = Protein.proteomeIdLookupTable.get(self.organismID)
    
    def __repr__(self):
        return f"Protein({self.accession}, {self.entryName}, {self.weight} kDa, AB={self.abundance})"
    
    def __str__(self):
        return f"Protein({self.accession}, {self.entryName}, {self.weight} kDa, AB={self.abundance})"
    
    # ------------------------------------------------------------------------------
    # ------------------------------ Internal methods ------------------------------
    def setSequenceDependentAttributes(self):
        """
        INTERNAL: recompute sequence-derived attributes after sequence changes.
        Ensure this method is called after any sequence mutation.
        :return: None. Updates weight, hydrophobicity, and isoelectric_point.

        """
        self.weight = self.calculate_weight(self.sequence)
        self.hydrophobicity = self.calculate_hydrophobicity(self.sequence)
        self.isoelectric_point = self.calculate_isoelectric_point(self.sequence)

    def calculate_weight(self,sequence):
        """
        Calculate the approximate molecular weight of a polypeptide sequence in kilodaltons (kDa).
        Parameters
        ----------
        sequence : str
            Amino acid sequence using single-letter codes. Residues not found in self.AMINO_ACID_WEIGHTS are treated as weight 0.
        Returns
        -------
        float
            Molecular weight in kDa, rounded to two decimal places. Computation sums residue masses from self.AMINO_ACID_WEIGHTS,
            subtracts 18.0153 Da per peptide bond (len(sequence) - 1) to account for released water molecules, and converts daltons to kDa.
        Notes
        -----
        - Empty or very short sequences will be handled (unknown residues contribute 0).
        """
        
        rawSequenceWeight = sum([self.AMINO_ACID_WEIGHTS.get(aa,0) for aa in sequence])
        waterWeight = 18.0153 * (len(sequence) - 1)
        sequenceWeight = rawSequenceWeight - waterWeight
        return round(sequenceWeight/1000, ndigits=N_DIGITS_GLOBAL)
    
    def calculate_hydrophobicity(self,sequence):
        """
        Calculate the average hydrophobicity of a polypeptide sequence using the Kyte-Doolittle scale.
        Parameters
        ----------
        sequence : str
            Amino acid sequence using single-letter codes. Residues not found in self.AMINO_ACID_HYDROPHOBICITY are treated as hydrophobicity 0.
        Returns
        -------
        float
            Average hydrophobicity, rounded to two decimal places. Computation sums residue hydrophobicities from self.AMINO_ACID_HYDROPHOBICITY
            and divides by sequence length.
        Notes
        -----
        - Empty sequences will return an average hydrophobicity of 0.0.
        """
        if len(sequence) == 0:
            return 0.0
        totalHydrophobicity = sum([self.AMINO_ACID_HYDROPHOBICITY.get(aa,0) for aa in sequence])
        averageHydrophobicity = totalHydrophobicity / len(sequence)
        return round(averageHydrophobicity, ndigits=N_DIGITS_GLOBAL)
    
    def calculate_isoelectric_point(self,sequence):
        """
        Calculate the isoelectric point (pI) of a polypeptide sequence.
        Parameters
        ----------
        sequence : str
            Amino acid sequence using single-letter codes.
        Returns
        -------
        float
            Isoelectric point (pI) of the sequence, rounded to two decimal places.
        Notes
        -----
        - Utilizes Biopython's SeqUtils module for pI calculation.
        - Empty sequences will return a pI of 0.0.
        """
        if len(sequence) == 0:
            return 0.0
        pI = IsoelectricPoint(str(sequence)).pi()
        return round(pI, ndigits=N_DIGITS_GLOBAL)
    
    
    def __cleaveSequence(self,startEnd,type="signal"):
        """
        INTERNAL: perform sequence cleavage (used by signal_peptide_removal module).

        Warning: class-specific; updates sequence, recalculates attributes, logs modification.

        :param startEnd: Tuple/list (start, end) 1-based indices to remove.
        :param type: Text label for modification log.
        :return: None.

        """
        if startEnd[0] and startEnd[1]:
            # -1 to convert to 0-based indexing
            startIndex = startEnd[0] - 1
            endIndex = startEnd[1] - 1
            # Remove sequence from startIndex to endIndex (inclusive)
            self.sequence = self.sequence[:startIndex] + self.sequence[endIndex + 1:]
            self.setSequenceDependentAttributes()
            self.modifications.append(f"Performed sequence cleavage of type: {type}")


    def __depleteAbundance(self,depletionPercentage):
        """
        INTERNAL: scale abundance by (1 - depletionPercentage).

        Warning: class-specific; used by affinity_depletion backend.

        :param depletionPercentage: Fraction between 0 and 1 to remove.
        :return: None. Mutates abundance and logs modification.

        """
        if self.abundance is not None:
            self.abundance *= (1 - depletionPercentage)
            self.abundance = max(self.abundance, 0)
            self.modifications.append(f"Depleted by {depletionPercentage*100}%")
            
    def get_fasta(self):
        """
        Return a FASTA-formatted string for this protein.
        Constructs a single-line FASTA header from instance getters in the form:
        ">DB|ACCESSION|ENTRYNAME proteinNamePlaceholder OS=Homo Sapiens OX=<organism_id> GN=<gene_name> PE=<protein_existence> SV=<sequence_version> AB=<abundance> MD=<True|False>"
        (MD is True if self.modifications is non-empty.)
        The protein sequence is wrapped to 60 characters per line.
            str: FASTA entry (header + newline-separated sequence lines).
        """
        header = f">{self.get_db()}|{self.get_accession()}|{self.get_entry_name()} \
        proteinNamePlaceholder OS=Homo Sapiens OX={self.get_organism_id()} \
        GN={self.get_gene_name()} PE={self.get_protein_existence()} \
        SV={self.get_sequence_version()} AB={self.get_abundance()} \
        MD={False if len(self.modifications) == 0 else True}"
        fasta_lines = [header]
        seq = str(self.sequence)
        for i in range(0, len(seq), 60):
            fasta_lines.append(seq[i:i+60])
        return "\n".join(fasta_lines)
    
    
    # ------------------------------- Set Functions ---------------------------------
    # -------------------------------------------------------------------------------

    def set_abundance(self, abundance):
        self.abundance = float(abundance)


    # ------------------------------- Get Functions ---------------------------------
    # -------------------------------------------------------------------------------
    
    def get_sequence(self):
        return self.sequence

    def get_weight(self):
        return getattr(self, "weight", None)
    
    def get_hydrophobicity(self):
        return getattr(self, "hydrophobicity", None)
    
    def get_isoelectricpoint(self):
        return getattr(self, "isoelectric_point", None)

    def get_db(self):
        return getattr(self, "db", None)

    def get_accession(self):
        return getattr(self, "accession", None)

    def get_entry_name(self):
        return getattr(self, "entryName", None)

    def get_protein_name(self):
        return getattr(self, "proteinName", None)

    def get_organism(self):
        return getattr(self, "organism", None)

    def get_organism_id(self):
        return getattr(self, "organismID", None)

    def get_gene_name(self):
        return getattr(self, "geneName", None)

    def get_protein_existence(self):
        return getattr(self, "proteinExistence", None)

    def get_sequence_version(self):
        return getattr(self, "sequenceVersion", None)

    def get_abundance(self):
        return getattr(self, "abundance", 0.0)

    def get_sequence_length(self):
        return len(self.sequence) if self.sequence is not None else 0
    
    # -------------------------------------------------------------------------------
    # ------------------------------- Get Functions ---------------------------------

    
    
def main():
    albuminHeader = ">sp|P02768|ALBU_HUMAN Albumin OS=Homo sapiens OX=9606 GN=ALB PE=1 SV=2 AB=40.0"
    albuminSequence = "MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLFFAKRYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAVARLSQRFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADDRADLAKYICENQDSISSKLKECCEKPLLEKSHCIAEVENDEMPADLPSLAADFVESKDVCKNYAEAKDVFLGMFLYEYARRHPDYSVVLLLRLAKTYETTLEKCCAAADPHECYAKVFDEFKPLVEEPQNLIKQNCELFEQLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPCAEDYLSVVLNQLCVLHEKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHADICTLSEKERQIKKQTALVELVKHKPKATKEQLKAVMDDFAAFVEKCCKADDKETCFAEEGKKLVAASQAALGL"
    tempObj = Protein(header=albuminHeader, sequence=albuminSequence)
    print(tempObj.get_weight())
    
if __name__ == '__main__':
    import time
    start_time = time.time()
    main()
    print('--- %s seconds ---' % (time.time() - start_time))