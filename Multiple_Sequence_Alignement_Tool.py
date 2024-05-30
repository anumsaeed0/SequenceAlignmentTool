import sys
import os
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QFrame, QProgressDialog, QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QTextEdit, QLabel, QLineEdit, QRadioButton, QGroupBox, QGridLayout, QMessageBox, QFileDialog
from PyQt5.QtGui import QFont
from Bio.Align import PairwiseAligner
from PyQt5 import QtCore
from PyQt5.QtWidgets import QMessageBox
from PyQt5.QtWidgets import QDialog
from PyQt5.QtSvg import QSvgWidget

from Bio import Phylo
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Align import MultipleSeqAlignment
from io import StringIO

from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import subprocess
import requests
import time
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

from pymsa import MSA, Entropy, PercentageOfNonGaps, PercentageOfTotallyConservedColumns, Star, SumOfPairs
from pymsa import PAM250, Blosum62, FileMatrix
from pymsa.util.fasta import print_alignment

from io import StringIO
from Bio import Phylo

class SeqAlignApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.initUI()
        self.child_windows = []  # List to track child windows

    def initUI(self):
        self.setWindowTitle("Multiple Sequence Alignment Tool")
        self.setGeometry(100, 100, 800, 600)

        # Main widget and layout
        widget = QWidget(self)
        self.setCentralWidget(widget)
        vbox = QVBoxLayout()

        # Heading
        heading_label = QLabel("Multiple Sequence Alignment Tool")
        heading_label.setStyleSheet("font-size: 22pt; font-weight: bold")
        heading_label.setAlignment(QtCore.Qt.AlignCenter)
        vbox.addWidget(heading_label)

        # Line
        line = QFrame()
        line.setFrameShape(QFrame.HLine)
        line.setFrameShadow(QFrame.Sunken)
        vbox.addWidget(line)

        # Sequence type selection
        type_groupbox = QGroupBox("Select Sequence Type")
        type_groupbox.setStyleSheet("font-size: 10pt;")
        type_layout = QHBoxLayout()
        self.dna_radio = QRadioButton("DNA")
        self.protein_radio = QRadioButton("Protein")
        self.dna_radio.setChecked(True)
        type_layout.addWidget(self.dna_radio)
        type_layout.addWidget(self.protein_radio)
        type_groupbox.setLayout(type_layout)
        vbox.addWidget(type_groupbox)

        # Sequence input
        seq1_label = QLabel("Sequences :")
        seq1_label.setStyleSheet("font-weight: bold; font-size: 10pt;")

        self.seq1_edit = QTextEdit()
        self.seq1_edit.setFixedHeight(130)
        self.seq1_edit.setPlaceholderText("Paste or upload sequences in FASTA format")

        # Upload buttons
        upload_button1 = QPushButton("Load Sequences")
        upload_button1.setStyleSheet("font-size: 10pt; color: white; background-color: #e88532; border: 1px solid #e88532; border-radius: 5px")
        upload_button1.setFixedHeight(130)
        upload_button1.setMinimumWidth(100)  # Set the minimum width here
        upload_button1.clicked.connect(lambda: self.upload_sequence(self.seq1_edit))

        # further styling
        seq_layout = QGridLayout()
        seq_layout.addWidget(seq1_label, 0, 0)
        seq_layout.addWidget(self.seq1_edit, 0, 1)
        seq_layout.addWidget(upload_button1, 0, 2)  # Upload button for Sequences

        vbox.addLayout(seq_layout)

        # Phylogenetic selection
        Tree_groupbox = QGroupBox("Create Phylogenetic Tree")
        Tree_groupbox.setStyleSheet("font-size: 10pt;")
        Tree_layout = QHBoxLayout()
        self.tree_yes_radio = QRadioButton("Yes")
        self.tree_no_radio = QRadioButton("No")
        self.tree_yes_radio.setChecked(True)

        Tree_layout.addWidget(self.tree_yes_radio)
        Tree_layout.addWidget(self.tree_no_radio)
        Tree_groupbox.setLayout(Tree_layout)
        vbox.addWidget(Tree_groupbox)

        # Alignment button
        align_button = QPushButton("Align Sequences")
        align_button.setStyleSheet("font-size: 12pt; color: white; background-color: #007bff; border: 1px solid #007bff; border-radius: 5px")
        align_button.setFixedHeight(50)
        align_button.clicked.connect(self.perform_alignment)

        reset_button = QPushButton("Reset")
        reset_button.setStyleSheet("font-size: 12pt; color: white; background-color: #dc3545; border: 1px solid #dc3545; border-radius: 5px")
        reset_button.setFixedHeight(50)
        reset_button.clicked.connect(self.reset_fields)

        hbox_buttons = QHBoxLayout()
        hbox_buttons.addWidget(align_button)
        hbox_buttons.addWidget(reset_button)
        vbox.addLayout(hbox_buttons)

        widget.setLayout(vbox)

    def perform_alignment(self):
        seq1 = self.seq1_edit.toPlainText().strip().upper()

        if seq1 == "":
            QMessageBox.critical(self, "Error", "Sequences (at least 3) are required!")
            return
        else:
            sequences = self.extractSeqs(seq1)
            if not self.validateSeq(sequences):
                QMessageBox.critical(self, "Error", "Invalid Sequences entered, Please enter a valid sequence!")
                return

        email = "anumsd01@gmail.com"

        # Displaying the message box while the job is running
        progress_dialog = QProgressDialog("Alignment job is running. Please wait...", None, 0, 0, self)
        progress_dialog.setWindowTitle("Alignment in progress")
        progress_dialog.setWindowModality(QtCore.Qt.WindowModal)
        progress_dialog.show()

        # Perform global alignment and generate phylogenetic tree
        alignment = self.run_clustal_omega(sequences, email)

        # Close the message box after the job is completed
        progress_dialog.close()

        if alignment:
            alignments = self.parse_alignment(alignment)
            alignments_seq = [(key, value) for key, value in alignments.items()]
            score_text, alignment_result = self.run_all_scores(alignments_seq) 
            
            tree = None
            if self.tree_yes_radio.isChecked():
                self.generate_phylogenetic_tree(alignments)

            # Display result in the new window
            app = QApplication(sys.argv)
            dialog = ResultDialog(score_text, alignment, self)
            dialog.show()
            dialog.exec_()
            sys.exit(app.exec_())
        else:
            QMessageBox.critical(self, "Error", "Alignment failed!")

    def upload_sequence(self, text_edit):
        file_dialog = QFileDialog()
        file_dialog.setNameFilter("Sequence files (*.txt *.fasta)")
        file_dialog.setViewMode(QFileDialog.List)
        file_dialog.setFileMode(QFileDialog.ExistingFile)

        if file_dialog.exec_():
            file_path = file_dialog.selectedFiles()[0]
            with open(file_path, 'r') as file:
                sequence = file.read()
                text_edit.setPlainText(sequence)
    
    def run_clustal_omega(self, sequences, email, datatype="protein"):
        # Set the endpoint URL for Clustal Omega
        url = "https://www.ebi.ac.uk/Tools/services/rest/clustalo/run/"

        # Define payload including sequences and email
        payload = {
            "email": email,
            "sequence": "\n".join([f">{name}\n{seq}" for name, seq in sequences.items()]),
            "outfmt": "clustal_num",
            "datatype": datatype,  # Change datatype here
            "action": "doAlignment"
        }

        # Send request to Clustal Omega web service
        response = requests.post(url, data=payload)

        # Check if request was successful
        if response.status_code == 200:
            # Retrieve job ID
            job_id = response.text.strip()
            print("Job ID:", job_id)

            # Wait for the job to complete
            while True:
                status_url = f"https://www.ebi.ac.uk/Tools/services/rest/clustalo/status/{job_id}"
                status_response = requests.get(status_url)
                status = status_response.text.strip()

                if status == "FINISHED":
                    break
                elif status == "ERROR":
                    print("Job encountered an error.")
                    return None
                else:
                    print("Job is still running...")
                    time.sleep(20) # Wait for 20 seconds before checking again
                            # Retrieve alignment results
            result_url = f"https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/{job_id}/aln-clustal_num"
            result_response = requests.get(result_url)
            alignment = result_response.text


            # Write alignment to a FASTA file
            alignment_file = f"alignment_{job_id}.fasta"
            with open(alignment_file, "w") as f:
                f.write(alignment)
            
            return alignment
            
        else:
            print("Error:", response.status_code, response.text)
            return None

    def parse_alignment(self, alignment_text):
        alignment_list = alignment_text.split("\n")

        alignment_list.pop(0)

        align_dict = dict()

        for alignment in alignment_list:
            if alignment.strip() != "":
                align_a = alignment.split()
                if len(align_a)==3:
                    seq_id = align_a[0]
                    seq_align = align_a[1]
                    if seq_id not in align_dict.keys():
                        align_dict[seq_id] = seq_align
                    else:
                        align_dict[seq_id]+=seq_align
                    align_dict[seq_id].replace(" ", "")

        return align_dict

    def extractSeqs(self, sequence):
        if sequence.startswith(">"):
            sequences = sequence.split("\n")

            MSeq = dict()
            id_seq = ""
            for seq in sequences:
                if seq.startswith(">"):
                    if seq[1:] not in MSeq.keys():
                        MSeq[seq[1:]] = ""
                        id_seq = seq[1:]
                else:
                    MSeq[id_seq]+=seq

            if len(MSeq) >= 3:
                return MSeq
            else:
                return None
        else:
            return None

        
    # Define the validateSeq method
    def validateSeq(self, sequences):
        if sequences is not None:
            sequence_type = "DNA" if self.dna_radio.isChecked() else "Protein"

            # Define valid characters based on the sequence type
            if sequence_type == 'DNA':
                valid_chars = set("ATGC")
            elif sequence_type == 'Protein':
                valid_chars = set("ACDEFGHIKLMNPQRSTVWY")
            
            # Loop through each sequence
            for seq in sequences.values():
                seq = seq.strip().upper()  # Convert sequence to uppercase and remove leading/trailing whitespace

                # Check if the sequence is empty
                if seq == "":
                    return False

                # Check if the sequence contains only valid characters
                if not set(seq).issubset(valid_chars):
                    return False

            # All sequences passed validation
            return True

    def run_all_scores(self, sequences: list) -> str:
        aligned_sequences = [pair[1] for pair in sequences]
        sequences_id = [pair[0] for pair in sequences]

        msa = MSA(aligned_sequences, sequences_id)
        
        score_text = ""

        # Percentage of non-gaps and totally conserved columns
        non_gaps = PercentageOfNonGaps(msa)
        totally_conserved_columns = PercentageOfTotallyConservedColumns(msa)

        percentage = non_gaps.compute()
        score_text += "Percentage of non-gaps: {0} %\n".format(percentage)

        conserved = totally_conserved_columns.compute()
        score_text += "Percentage of totally conserved columns: {0}\n".format(conserved)

        # Entropy
        value = Entropy(msa).compute()
        score_text += "Entropy score: {0}\n".format(value)

        # Sum of pairs
        value = SumOfPairs(msa, Blosum62()).compute()
        score_text += "Sum of Pairs score (Blosum62): {0}\n".format(value)

        value = SumOfPairs(msa, PAM250()).compute()
        score_text += "Sum of Pairs score (PAM250): {0}\n".format(value)

        # Get the directory where the main script is located
        script_dir = os.path.dirname(os.path.abspath(__file__))

        # Construct the absolute path for the PAM380.txt file
        pam380_file_path = os.path.join(script_dir, "PAM380.txt")

        # Use the constructed file path in the FileMatrix initialization
        value = SumOfPairs(msa, FileMatrix(pam380_file_path)).compute()
        score_text += "Sum of Pairs score (PAM380): {0}\n".format(value)

        # Star
        value = Star(msa, Blosum62()).compute()
        score_text += "Star score (Blosum62): {0}\n".format(value)

        value = Star(msa, PAM250()).compute()
        score_text += "Star score (PAM250): {0}\n".format(value)

        return score_text, msa
    def generate_phylogenetic_tree(self, alignment_data):
        # Parse the alignment data
        alignment = MultipleSeqAlignment([])
        for key, value in alignment_data.items():
            if value.strip():
                seq_record = SeqRecord(Seq(value), id=key)
                alignment.append(seq_record)

        # Calculate the distance matrix
        calculator = DistanceCalculator('identity')
        distance_matrix = calculator.get_distance(alignment)

        # Build the tree
        constructor = DistanceTreeConstructor()
        tree = constructor.upgma(distance_matrix)

        # Draw and display the tree
        Phylo.draw(tree)


    def reset_fields(self):
        self.seq1_edit.clear()
        # Reset radio buttons to default
        self.dna_radio.setChecked(True)  # Resets to DNA being selected
        self.tree_yes_radio.setChecked(True) 

        for window in self.child_windows:
            window.close()
        self.child_windows.clear()


class ResultDialog(QDialog):
    def __init__(self, score, alignment_result, parent=None):
        super().__init__(parent)
        self.initUI(score, alignment_result)

    def initUI(self, score, alignment_result):
        self.setWindowTitle("Alignment Result")
        self.setGeometry(200, 200, 800, 600)

        vbox = QVBoxLayout()

        score_label = QLabel("Alignment Score:")
        vbox.addWidget(score_label)

        self.score_edit = QTextEdit()
        self.score_edit.setReadOnly(True)
        self.score_edit.setStyleSheet("font-size: 12pt; font:'Courier New';")  # Set the font size here
        self.score_edit.setPlainText(score)
        vbox.addWidget(self.score_edit)

        alignment_label = QLabel("Alignment Result:")
        vbox.addWidget(alignment_label)

        self.alignment_edit = QTextEdit()
        font = QFont("Courier New")
        font.setStyleHint(QFont.Monospace)  # Optional: hint to use monospaced style
        self.alignment_edit.setReadOnly(True)
        self.alignment_edit.setStyleSheet("font-size: 12pt;")  # Set the font size here
        self.alignment_edit.setFont(font)
        self.alignment_edit.setPlainText(alignment_result)
        vbox.addWidget(self.alignment_edit)
        
        # Save button
        self.save_button = QPushButton("Save Alignment Result")
        self.save_button.setStyleSheet("font-size: 12pt; color: white; background-color: #28a745; border: 1px solid #28a745; border-radius: 5px")
        self.save_button.setFixedHeight(30)
        self.save_button.clicked.connect(self.save_results)
        vbox.addWidget(self.save_button)
        
        # Set window flags
        self.setWindowFlags(Qt.Window | Qt.WindowMinimizeButtonHint | Qt.WindowMaximizeButtonHint | Qt.WindowCloseButtonHint)
        
        self.setLayout(vbox)


    def save_results(self):
        file_path, _ = QFileDialog.getSaveFileName(self, "Save File", "", "Text Files (*.txt);;All Files (*)")
        if file_path:
            with open(file_path, 'w') as file:
                file.write(self.score_edit.toPlainText())
                file.write(self.alignment_edit.toPlainText())
                QMessageBox.information(self, "Save Successful", "The alignment results were successfully saved!")



def main():
    app = QApplication(sys.argv)
    ex = SeqAlignApp()
    ex.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()