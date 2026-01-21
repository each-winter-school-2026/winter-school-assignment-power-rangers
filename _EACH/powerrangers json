import json
from SDS.SDS_PAGE import getExampleSDS_PAGE,virtualSDSPage_2DGaussian,parseFasta
from _EACH.protein import Protein
from utils.helperFunctions import extractSetting

def select(moduleIdentifier,selectedSettings,moduleData):
    """
    Dispatch the requested module to its backend implementation and return the simulated SDS-PAGE image.

    Settings (see module JSON):
    - Varies per module; each module function documents its own settings.

    :param moduleIdentifier: String module id (e.g., "fasta_input") coming from the frontend JSON definition.
    :param selectedSettings: Dict of user-selected values for the active module.
    :param moduleData: Loaded JSON defining all modules and their settings (options, defaults, etc.).
    :return: SDS-PAGE plot generated from the proteins produced by the module.

    """
    match moduleIdentifier:
        case "fasta_input":
            proteins = fasta_input(moduleIdentifier,selectedSettings,moduleData)
            return virtualSDSPage_2DGaussian(proteins)
        case "affinity_depletion":
            proteins = affinity_depletion(moduleIdentifier,selectedSettings,moduleData)
            return virtualSDSPage_2DGaussian(proteins)
        case "molecular_weight_cutoff":
            proteins = molecular_weight_cutoff(moduleIdentifier,selectedSettings,moduleData)
            return virtualSDSPage_2DGaussian(proteins)
        case "signal_peptide_removal":
            proteins = signal_peptide_removal(moduleIdentifier,selectedSettings,moduleData)
            return virtualSDSPage_2DGaussian(proteins)
        case "isoelectric_focussing":
            proteins = isoelectric_focussing(moduleIdentifier,selectedSettings,moduleData)
            return virtualSDSPage_2DGaussian(proteins)
        case "reversed_phase_chromatography":
            proteins = reversed_phase_chromatography(moduleIdentifier,selectedSettings,moduleData)
            return virtualSDSPage_2DGaussian(proteins)
        case "SDS_page_fractionation":
            proteins = SDS_page_fractionation(moduleIdentifier,selectedSettings,moduleData)
            return virtualSDSPage_2DGaussian(proteins)
        case "unique_module_identifier":
            proteins = newModule(moduleIdentifier,selectedSettings,moduleData)
            return virtualSDSPage_2DGaussian(proteins)
        case "example_module":
            # Does not perform any real processing; for demonstration only.
            proteins = exampleModule(moduleIdentifier,selectedSettings,moduleData)
            return virtualSDSPage_2DGaussian(proteins)
        case _: # Add new modules above 
            # Do not add modules below
            raise NotImplementedError(f"Module: {moduleIdentifier} is not implemented yet.")
        

def fasta_input(moduleIdentifier, selectedSettings,moduleData):
    """
    Load proteins from a selected FASTA file and convert sequences into Protein objects.

    Settings (fasta_input.json):
    - "Select FASTA file" (ChoiceField): maps a human-readable name to a FASTA filepath.

    :param selectedSettings: Dict containing the chosen FASTA file key from the UI.
    :param moduleData: Module definitions used to map the chosen key to an actual FASTA path.
    :return: List of Protein instances created from the FASTA sequences.

    """
    # Resolve selected label to actual FASTA path
    filePath = extractSetting("Select FASTA file",moduleIdentifier,selectedSettings,moduleData)
    sequences = parseFasta(filePath)
    Protein.deleteAllProteins()
    proteinList = []
    for header, sequence in sequences.items():
        proteinList.append(Protein(header, sequence))
    return proteinList

def affinity_depletion(moduleIdentifier, selectedSettings,moduleData):
    """
    Apply immunoaffinity depletion to targeted genes or gene groups, optionally with a custom percentage.

    Settings (affinity_depletion.json):
    - "Depletion Target" (ChoiceField): preset kits mapping to gene lists with default depletion %.
    - "Use custom depletion percentage" (BooleanField): toggles custom % for all selected targets.
    - "Custom Depletion Percentage" (DecimalField): custom percentage (0-100) applied when enabled.

    :param moduleIdentifier: Identifier for the current module.
    :param selectedSettings: Dict with depletion target selection and optional custom depletion percentage.
    :param moduleData: Module definitions providing target options and gene-group mappings.
    :return: Updated list of Protein objects after depletion.

    """
    # ChoiceField: resolve kit name to list of (target, default%) tuples
    depletionTargets = extractSetting("Depletion Target",moduleIdentifier,selectedSettings,moduleData)
    overwriteDepletionPercentage = extractSetting("Use custom depletion percentage",moduleIdentifier,selectedSettings,moduleData)
    customDepletionPercentage = extractSetting("Custom Depletion Percentage",moduleIdentifier,selectedSettings,moduleData)/100.0
    geneGroups = json.load(open('modules/geneGroups.json'))
    genesToDeplete = {}
    
    for target,depletionPercentage in depletionTargets:
        print(target,depletionPercentage)
        if target.startswith("Group:"):
            groupName = target.split("Group:")[1]
            # Expand group into individual genes
            for gene in geneGroups[groupName]:
                genesToDeplete[gene] = customDepletionPercentage if overwriteDepletionPercentage else depletionPercentage
        else:
            genesToDeplete[target] = customDepletionPercentage if overwriteDepletionPercentage else depletionPercentage

    Protein.proteinImmunoaffinityDepletion(genesToDeplete)
    return Protein.getAllProteins()

def molecular_weight_cutoff(moduleIdentifier, selectedSettings,moduleData):
    """
    Filter proteins by molecular weight, keeping those above or below a user-defined cutoff.

    Settings (molecular_weight_cutoff.json):
    - "Weight Cutoff (kDa)" (DecimalField): numeric threshold.
    - "Keep Below/Above Cutoff" (ChoiceField): maps UI labels to "below" or "above" behavior.

    :param moduleIdentifier: Identifier for the current module.
    :param selectedSettings: Dict with cutoff value and whether to keep proteins above or below it.
    :param moduleData: Module definitions supplying the mapping for the keep-above/below choice.
    :return: Updated list of Protein objects after weight-based filtering.

    """
    weight_cutoff = extractSetting("Weight Cutoff (kDa)",moduleIdentifier,selectedSettings,moduleData)
    # ChoiceField: resolve label to behavior string
    keep_option = extractSetting("Keep Below/Above Cutoff",moduleIdentifier,selectedSettings,moduleData)
    if keep_option == "above":
        Protein.depleteProteinsByWeight(minWeight=weight_cutoff,maxWeight=None)
    elif keep_option == "below":
        Protein.depleteProteinsByWeight(minWeight=None,maxWeight=weight_cutoff)
    else: 
        raise ValueError(f"Invalid option for Keep Below/Above Cutoff: {keep_option}")
    return Protein.getAllProteins()


def signal_peptide_removal(moduleIdentifier, selectedSettings,moduleData):
    """
    Remove signal peptides from proteins according to the selected database configuration.

    Settings (signal_peptide_removal.json):
    - "Database" (ChoiceField): selects the database key (e.g., "uniprot") used by the cleavage logic.

    :param moduleIdentifier: Identifier for the current module.
    :param selectedSettings: Dict indicating which database option to use for cleavage.
    :param moduleData: Module definitions used to resolve the chosen database option.
    :return: Updated list of Protein objects after signal peptide cleavage.

    """
    # ChoiceField: resolve label to database key
    database = extractSetting("Database",moduleIdentifier,selectedSettings,moduleData)
    Protein.signalPeptideCleavage()
    return Protein.getAllProteins()


def isoelectric_focussing(moduleIdentifier, selectedSettings,moduleData):
    """
    Fractionate proteins by isoelectric point, keeping those inside or outside a specified pI range.

    Settings (isoelectric_focussing.json):
    - "Keep inside/outside isoelectric point range" (ChoiceField): maps to "inside" or "outside" behavior.
    - "Minimum pI" (DecimalField): lower bound of pI window.
    - "Maximum pI" (DecimalField): upper bound of pI window.

    :param moduleIdentifier: Identifier for the current module.
    :param selectedSettings: Dict with min/max pI values and the keep-inside/outside selection.
    :param moduleData: Module definitions resolving the keep-inside/outside option mapping.
    :return: Updated list of Protein objects after pI-based fractionation.

    """
    keepInsideOutside = extractSetting("Keep inside/outside isoelectric point range",moduleIdentifier,selectedSettings,moduleData)
    pI_min = extractSetting("Minimum pI",moduleIdentifier,selectedSettings,moduleData)
    pI_max = extractSetting("Maximum pI",moduleIdentifier,selectedSettings,moduleData)
    Protein.fractionateProteinsByIsoelectricPoint(keepInsideOutsideSelection=keepInsideOutside,minPI=pI_min,maxPI=pI_max)
    return Protein.getAllProteins()
    
def reversed_phase_chromatography(moduleIdentifier, selectedSettings,moduleData):
    """
    Fractionate proteins by hydrophobicity, keeping those inside or outside a specified hydrophobicity range.

    Settings (reversed_phase_chromatography.json):
    - "Keep inside/outside hydrophobicity range" (ChoiceField): maps to "inside" or "outside" behavior.
    - "Minimum hydrophobicity" (DecimalField): lower bound of hydrophobicity window.
    - "Maximum hydrophobicity" (DecimalField): upper bound of hydrophobicity window.

    :param moduleIdentifier: Identifier for the current module.
    :param selectedSettings: Dict with min/max hydrophobicity values and the keep-inside/outside selection.
    :param moduleData: Module definitions resolving the keep-inside/outside option mapping.
    :return: Updated list of Protein objects after hydrophobicity-based fractionation.

    """
    keepInsideOutside = extractSetting("Keep inside/outside hydrophobicity range",moduleIdentifier,selectedSettings,moduleData)
    hydrophobicity_min = extractSetting("Minimum hydrophobicity",moduleIdentifier,selectedSettings,moduleData)
    hydrophobicity_max = extractSetting("Maximum hydrophobicity",moduleIdentifier,selectedSettings,moduleData)
    Protein.fractionateProteinsByHydrophobicity(keepInsideOutsideSelection=keepInsideOutside,minHydrophobicity=hydrophobicity_min,maxHydrophobicity=hydrophobicity_max)
    return Protein.getAllProteins()
    
def SDS_page_fractionation(moduleIdentifier, selectedSettings,moduleData):
    """
    Fractionate proteins by molecular weight, keeping those inside or outside a specified weight range.

    Settings (SDS_page_fractionation.json):
    - "Keep inside/outside weight range" (ChoiceField): maps to "inside" or "outside" behavior.
    - "Minimum weight (kDa)" (DecimalField): lower bound of weight window.
    - "Maximum weight (kDa)" (DecimalField): upper bound of weight window.

    :param moduleIdentifier: Identifier for the current module.
    :param selectedSettings: Dict with min/max weight values and the keep-inside/outside selection.
    :param moduleData: Module definitions resolving the keep-inside/outside option mapping.
    :return: Updated list of Protein objects after weight-based fractionation.

    """
    # ChoiceField: resolve label to behavior string
    keepInsideOutside = extractSetting("Keep inside/outside molecular weight range",moduleIdentifier,selectedSettings,moduleData)
    weight_min = extractSetting("Minimum weight (kDa)",moduleIdentifier,selectedSettings,moduleData)
    weight_max = extractSetting("Maximum weight (kDa)",moduleIdentifier,selectedSettings,moduleData)
    Protein.fractionateProteinsByMolecularWeight(keepInsideOutsideSelection=keepInsideOutside,minWeight=weight_min,maxWeight=weight_max)
    return Protein.getAllProteins()
    
    
def newModule(moduleIdentifier,selectedSettings,moduleData):
    # The first step is to access the settings chosen by the user. 
    
    # We start by extracting the ChoiceField. This will be resolved to the internal values within the json file. 
    choiceFieldChosenOption = extractSetting("A field where the user can choose from multiple options",moduleIdentifier,selectedSettings,moduleData)
    # Next we extract the DecimalField where the user can enter their own number.
    chosenNumber = extractSetting("A field where the user can enter their own number",moduleIdentifier,selectedSettings,moduleData)
    # Next lets loop through all proteins and set their abundances to the user supplied values. 
    
    # First get all the proteins in a nice list
    proteins = Protein.getAllProteins()

    # Now lets loop through all the proteins
    for protein in proteins:
        # And set their abundance to the chosen number
        protein.set_abundance(chosenNumber)
    return proteins

def exampleModule(moduleIdentifier,selectedSettings,moduleData):
    # This is an example module that does nothing. It simply returns all proteins as is.
    singleChoiceFieldOption = extractSetting("Single choice field",moduleIdentifier,selectedSettings,moduleData)
    multipleChoiceFieldOptions = extractSetting("Multiple choice field",moduleIdentifier,selectedSettings,moduleData)
    decimalFieldValue = extractSetting("Decimal field",moduleIdentifier,selectedSettings,moduleData)
    booleanFieldValue = extractSetting("Boolean field",moduleIdentifier,selectedSettings,moduleData)
    charFieldValue = extractSetting("Character field",moduleIdentifier,selectedSettings,moduleData)
    
    print(f"Single choice field selected option: {singleChoiceFieldOption}")
    print(f"Multiple choice field selected options: {multipleChoiceFieldOptions}")
    print(f"Decimal field value: {decimalFieldValue}")
    print(f"Boolean field value: {booleanFieldValue}")
    print(f"Character field value: {charFieldValue}")
    
    
    return Protein.getAllProteins()


