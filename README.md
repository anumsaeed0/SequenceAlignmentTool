# Sequence Alignment Tool

## Introduction
This tool provides a graphical user interface (GUI) for performing both pairwise and multiple sequence alignment (MSA) of biological sequences, such as DNA, RNA, and protein sequences. It is developed using Python and PyQt5 for the GUI components, and it utilizes Biopython for sequence alignment algorithms and Bioinformatics tools.

## Features
### Pairwise Sequence Alignment
- **Supported Sequences**: Supports alignment of DNA, RNA, and protein sequences.
- **Alignment Types**: Provides options for both global and local pairwise sequence alignment.
- **Customizable Penalties**: Allows users to specify match, mismatch, gap open, and gap extend penalties.
- **Sequence Upload**: Enables users to upload sequences from text or FASTA files.
- **Alignment Display**: Displays alignment results with customizable font and HTML formatting.
- **Result Saving**: Allows users to save alignment results to a text file.

### Multiple Sequence Alignment (MSA)
- **Sequence Input**: Accepts multiple sequences in FASTA format for MSA.
- **Alignment Tool Integration**: Integrates with Clustal Omega for performing MSA via web service.
- **Phylogenetic Tree Generation**: Generates a phylogenetic tree based on the MSA results.
- **Alignment Quality Metrics**: Calculates and displays various alignment quality metrics, including percentage of non-gaps, percentage of totally conserved columns, entropy score, sum of pairs scores, and star scores.
- **Customizable Score Matrix**: Allows users to choose from pre-defined scoring matrices (e.g., PAM250, BLOSUM62) or specify a custom matrix file.

## Installation
1. **Clone Repository**: Clone the repository to your local machine:
    ```bash
    git clone https://github.com/anumsaeed0/alignment-tool.git
    ```

2. **Install Dependencies**: Install the required Python packages using pip:
    ```bash
    pip install PyQt5
    pip install biopython
    pip install pymsa
    pip install requests
    pip install matplotlib
    ```

3. **Run the Application**: Execute the main Python script to launch the application:
    ```bash
    python Main_Window_Tool.py
    ```

## Usage
1. **Pairwise Sequence Alignment**:
    - Launch the application.
    - Select the type of sequences (DNA, RNA, or Protein).
    - Input or upload sequences for alignment.
    - Choose the alignment type (Global or Local).
    - Specify alignment penalties if needed.
    - Click the "Align Sequences" button to perform alignment.
    - View and save alignment results.

2. **Multiple Sequence Alignment (MSA)**:
    - Follow the steps for pairwise alignment to input or upload sequences.
    - Enable the "Create Phylogenetic Tree" option if required.
    - Click the "Align Sequences" button to perform MSA.
    - View alignment results, including the phylogenetic tree and alignment quality metrics.
    - Save alignment results if needed.

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments
- Biopython: For providing the tools for sequence alignment.
- PyQt5: For creating the graphical user interface.
- Clustal Omega: For the web service used in multiple sequence alignment.
