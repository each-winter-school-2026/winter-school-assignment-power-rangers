[![Open in Codespaces](https://classroom.github.com/assets/launch-codespace-2972f46106e565e64193e422d61a12cf1da4916b45550586e14ef0a7c637dd04.svg)](https://classroom.github.com/open-in-codespaces?assignment_repo_id=22297801)
# EWOKS

**Experimental Workflow Optimization through Knowledge Synthesis**

EWOKS is a Django-based web application for simulating protein purification and fractionation workflows. It combines interactive web interface design, computational modeling of protein properties, and virtual gel electrophoresis visualization to create an educational platform for proteomics learning and workflow prototyping.

---

## About This Repository

This is a sandbox environment for developing and testing EWOKS components:
- **Web Interface**: Django application with interactive module configuration
- **Backend System**: Modular processing pipeline with JSON-based module definitions
- **Protein Simulation**: Object-oriented protein class with physicochemical property calculations
- **Visualization**: Virtual SDS-PAGE gel rendering engine

The repository contains both the complete application code and comprehensive documentation for extending the system with new modules.

---

## Repository Structure

```
.
├── app/                              # Django web application
│   ├── EWOKS/                        # Django project configuration
│   │   ├── settings.py               # Django settings
│   │   ├── urls.py                   # URL routing
│   │   └── wsgi.py                   # WSGI configuration
│   ├── EWOKS_Interface/              # Main Django app
│   │   ├── views.py                  # Request handlers
│   │   ├── models.py                 # Database models
│   │   ├── formConstructor.py        # Form generation
│   │   ├── templates/                # HTML templates
│   │   └── static/                   # CSS, JavaScript, images
│   └── manage.py                     # Django management script
│
├── _EACH/                            # Backend module system (core)
│   ├── protein.py                    # Protein class definition
│   ├── modules.py                    # Module implementations
│   ├── protein_reference.md          # Protein class API documentation
│   ├── HOW_TO_IMPLEMENT_A_MODULE.md # Step-by-step module creation guide
│   └── modules/                      # Module JSON configurations
│       ├── _MODULE_TEMPLATE.md       # JSON field type reference
│       ├── fasta_input.json
│       ├── molecular_weight_cutoff.json
│       ├── affinity_depletion.json
│       ├── signal_peptide_removal.json
│       ├── isoelectric_focussing.json
│       ├── reversed_phase_chromatography.json
│       └── SDS_page_fractionation.json
│
├── SDS/                              # SDS-PAGE simulation engine
│   └── SDS_PAGE.py                   # Gel visualization and rendering
│
├── data/                             # Proteome data files
│   ├── proteomes/                    # Species proteome files
│   │   ├── Human/
│   │   └── Arabidopsis/
│   └── signalProteomes/              # Signal peptide databases
│
├── modules/                          # Module utilities
│   └── signal.py                     # Signal peptide lookup system
│
├── utils/                            # Utility functions
│   ├── helperFunctions.py            # Shared helpers (extractSetting, etc.)
│   ├── parseFasta.py                 # FASTA parsing
│   └── appendAbundance.py            # Abundance data processing
│
├── examples/                         # Example datasets
│   └── human_plasma_proteins*.tsv    # Plasma proteome data
│
├── requirements.txt                  # Python dependencies
└── README.md                         # This file
```

---

## Core Components

### Backend Module System (`_EACH/`)

The modular architecture allows researchers to build complex protein purification workflows by chaining together analysis steps. Each module:
- Defines user-facing parameters in JSON format
- Implements processing logic in Python
- Operates on the global Protein registry
- Returns modified protein population to the visualization system

**Available Modules:**
- **FASTA Input**: Load protein sequences with abundance data
- **Affinity Depletion**: Remove proteins by name with configurable depletion efficiency
- **Molecular Weight Cutoff**: Filter proteins above/below specified weight thresholds
- **Signal Peptide Removal**: Cleave N-terminal signal sequences
- **Isoelectric Focusing**: Separate proteins by pI (isoelectric point)
- **Reversed Phase Chromatography**: Fractionate by hydrophobicity
- **SDS-PAGE Fractionation**: Simulate gel-based molecular weight separation

### Protein Class (`_EACH/protein.py`)

Represents individual proteins with:
- **Sequence data**: Bio.Seq object with full amino acid sequence
- **Physicochemical properties**: Calculated on instantiation
  - Molecular weight (with water loss deduction)
  - Isoelectric point (Biopython)
  - Hydrophobicity (Kyte-Doolittle scale)
- **Metadata**: UniProt fields (accession, gene name, organism, etc.)
- **Abundance**: Relative abundance value for SDS-PAGE visualization
- **Modification tracking**: Log of applied transformations

All proteins are automatically registered in a global class-level dictionary for easy access across modules.

### Visualization Engine (`SDS/SDS_PAGE.py`)

Generates simulated 2D gel images showing protein separation based on molecular weight and relative abundance. Used by the web interface to visualize workflow results.

### Web Interface (`app/EWOKS_Interface/`)

Django application providing:
- Interactive module configuration forms
- Workflow builder UI
- Real-time SDS-PAGE gel visualization
- Database persistence for workflows

---

## Key Technologies

| Component | Technology | Version |
|-----------|-----------|---------|
| Web Framework | Django | 5.2.8 |
| Biology | Biopython | 1.86 |
| Numerics | NumPy | 2.3.4 |
| Visualization | Matplotlib | 3.10.7 |
| Forms | django-widget-tweaks | - |

See [requirements.txt](requirements.txt) for complete dependency list.

---

## Documentation

### For Module Developers

Start here to learn how to create new analysis modules:
- [**HOW_TO_IMPLEMENT_A_MODULE.md**](_EACH/HOW_TO_IMPLEMENT_A_MODULE.md) - Complete step-by-step tutorial with practical example
- [**_MODULE_TEMPLATE.md**](_EACH/modules/_MODULE_TEMPLATE.md) - Reference for JSON field types and backend integration patterns
- [**protein_reference.md**](_EACH/protein_reference.md) - API documentation for the Protein class

### For Application Users

The web interface at `http://localhost:8000` provides an interactive guide to building workflows.

---

## Quick Start

### Installation

```bash
# Set PYTHONPATH to the top-level directory (required for imports)
# On Windows (PowerShell):
$env:PYTHONPATH = "C:\path\to\sandbox"

# On Windows (Command Prompt):
set PYTHONPATH=C:\path\to\sandbox

# On Linux/macOS:
export PYTHONPATH=/path/to/sandbox

# Install dependencies
pip install -r requirements.txt

# Start Django development server
python app/manage.py runserver
```

Access the application at `http://localhost:8000`

> **Note:** The PYTHONPATH variable allows Python to import modules from the `_EACH`, `SDS`, and `utils` directories. Adjust the path according to your installation directory.
### Building a Workflow

1. Load protein data using the FASTA Input module
2. Add analysis modules in sequence (e.g., Affinity Depletion → Molecular Weight Cutoff)
3. Configure each module's parameters through the web form
4. View simulated SDS-PAGE gel after each step
5. Export or save your workflow

---

## Data Sources

- **Human proteome**: UniProt curated data
- **Plasma abundance**: Human plasma proteome reference datasets
- **Signal peptides**: Pre-computed signal peptide predictions for UniProt proteins

---

## Development

This is an active sandbox environment for EWOKS development. Code structure and APIs may change as the project evolves.

For questions about implementation details, see the documentation in [`_EACH/`](_EACH/) directory.
