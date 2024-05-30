import sys
from PyQt5.QtWidgets import QFrame, QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QTextEdit, QLabel, QLineEdit, QRadioButton, QGroupBox, QGridLayout, QMessageBox, QFileDialog
from Bio.Align import PairwiseAligner
from PyQt5 import QtCore

class SeqAlignApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        self.setWindowTitle("Pairwise Sequence Alignment Tool")
        self.setGeometry(100, 100, 800, 600)

        # Main widget and layout
        widget = QWidget(self)
        self.setCentralWidget(widget)
        vbox = QVBoxLayout()

        # Heading
        heading_label = QLabel("Pairwise Sequence Alignment Tool")
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
        self.rna_radio = QRadioButton("RNA")
        self.protein_radio = QRadioButton("Protein")
        self.dna_radio.setChecked(True)
        type_layout.addWidget(self.dna_radio)
        type_layout.addWidget(self.rna_radio)
        type_layout.addWidget(self.protein_radio)
        type_groupbox.setLayout(type_layout)
        vbox.addWidget(type_groupbox)

        # Sequence input
        seq1_label = QLabel("Sequence 1:")
        seq1_label.setStyleSheet("font-weight: bold; font-size: 10pt;")
        seq2_label = QLabel("Sequence 2:")
        seq2_label.setStyleSheet("font-weight: bold; font-size: 10pt;")
        self.seq1_edit = QTextEdit()
        self.seq2_edit = QTextEdit()

        # Set fixed height for sequence input boxes
        self.seq1_edit.setFixedHeight(60)
        self.seq2_edit.setFixedHeight(60)

        # set placeholder text
        self.seq1_edit.setPlaceholderText("Paste or upload sequence")
        self.seq2_edit.setPlaceholderText("Paste or upload sequence")

        # Upload buttons
        upload_button1 = QPushButton("Load Sequence")
        upload_button1.setStyleSheet("font-size: 10pt; color: white; background-color: #e88532; border: 1px solid #e88532; border-radius: 5px")
        upload_button1.setFixedHeight(60)
        upload_button1.setMinimumWidth(120)  # Set the minimum width here
        upload_button1.clicked.connect(lambda: self.upload_sequence(self.seq1_edit))
        
        upload_button2 = QPushButton("Load Sequence")
        upload_button2.setStyleSheet("font-size: 10pt; color: white; background-color: #e88532; border: 1px solid #e88532; border-radius: 5px")
        upload_button2.setFixedHeight(60)
        upload_button2.setMinimumWidth(120)  # Set the minimum width here
        upload_button2.clicked.connect(lambda: self.upload_sequence(self.seq2_edit))

        # further styling
        seq_layout = QGridLayout()
        seq_layout.addWidget(seq1_label, 0, 0)
        seq_layout.addWidget(self.seq1_edit, 0, 1)
        seq_layout.addWidget(upload_button1, 0, 2)  # Upload button for Sequence 1
        
        seq_layout.addWidget(seq2_label, 1, 0)
        seq_layout.addWidget(self.seq2_edit, 1, 1)
        seq_layout.addWidget(upload_button2, 1, 2)  # Upload button for Sequence 2


        vbox.addLayout(seq_layout)

        # Alignment type selection
        align_groupbox = QGroupBox("Select Alignment Type")
        align_groupbox.setStyleSheet("font-size: 10pt;")
        align_layout = QHBoxLayout()
        self.global_radio = QRadioButton("Global Alignment")
        self.local_radio = QRadioButton("Local Alignment")
        self.global_radio.setChecked(True)
        align_layout.addWidget(self.global_radio)
        align_layout.addWidget(self.local_radio)
        align_groupbox.setLayout(align_layout)
        vbox.addWidget(align_groupbox)
        
        # Penalties input
        penalties_groupbox = QGroupBox("Penalties")
        penalties_groupbox.setStyleSheet("font-size: 10pt;")
        penalties_layout = QGridLayout()

        self.match_penalty_edit = QLineEdit("1")
        self.mismatch_penalty_edit = QLineEdit("-1")
        self.gap_open_penalty_edit = QLineEdit("-0.5")
        self.gap_extend_penalty_edit = QLineEdit("-0.1")

        match_label = QLabel("Match Penalty:")
        match_label.setStyleSheet("font-weight: bold;  font-size: 10pt;")
        mismatch_label = QLabel("Mismatch Penalty:")
        mismatch_label.setStyleSheet("font-weight: bold; font-size: 10pt;")
        gap_open_label = QLabel("Gap Open Penalty:")
        gap_open_label.setStyleSheet("font-weight: bold; font-size: 10pt;")
        gap_extend_label = QLabel("Gap Extend Penalty:")
        gap_extend_label.setStyleSheet("font-weight: bold; font-size: 10pt;")

        penalties_layout.addWidget(match_label, 0, 0)
        penalties_layout.addWidget(self.match_penalty_edit, 0, 1)
        penalties_layout.addWidget(mismatch_label, 0, 2)
        penalties_layout.addWidget(self.mismatch_penalty_edit, 0, 3)
        penalties_layout.addWidget(gap_open_label, 0, 4)
        penalties_layout.addWidget(self.gap_open_penalty_edit, 0, 5)
        penalties_layout.addWidget(gap_extend_label, 0, 6)
        penalties_layout.addWidget(self.gap_extend_penalty_edit, 0, 7)

        penalties_groupbox.setLayout(penalties_layout)
        vbox.addWidget(penalties_groupbox)

        # Alignment button
        align_button = QPushButton("Align Sequences")
        align_button.setStyleSheet("font-size: 12pt; color: white; background-color: #007bff; border: 1px solid #007bff; border-radius: 5px")
        align_button.setFixedHeight(30)
        align_button.clicked.connect(self.perform_alignment)

        reset_button = QPushButton("Reset")
        reset_button.setStyleSheet("font-size: 12pt; color: white; background-color: #dc3545; border: 1px solid #dc3545; border-radius: 5px")
        reset_button.setFixedHeight(30)
        reset_button.clicked.connect(self.reset_fields)

        hbox_buttons = QHBoxLayout()
        hbox_buttons.addWidget(align_button)
        hbox_buttons.addWidget(reset_button)
        vbox.addLayout(hbox_buttons)

        # Save button
        self.save_button = QPushButton("Save Alignment Result")
        self.save_button.setStyleSheet("font-size: 12pt; color: white; background-color: #28a745; border: 1px solid #28a745; border-radius: 5px")
        self.save_button.setFixedHeight(30)
        self.save_button.clicked.connect(self.save_results)
        self.save_button.hide()  # Hide initially
        hbox_buttons.addWidget(self.save_button)


        # Result display
        self.result_label = QLabel("Alignment Result:")
        self.result_label.setStyleSheet("font-weight: bold; font-size: 12pt;")
        self.result_label.hide()  # Hide the label initially
        vbox.addWidget(self.result_label)

        self.result_edit = QTextEdit()
        self.result_edit.setReadOnly(True)
        self.result_edit.setStyleSheet("font-size: 12pt;")  # Set the font size here
        self.result_edit.hide()  # Hide initially
        vbox.addWidget(self.result_edit)

        widget.setLayout(vbox)

    def perform_alignment(self):
        seq1 = self.seq1_edit.toPlainText().replace('\n', '').strip().upper()
        seq2 = self.seq2_edit.toPlainText().replace('\n', '').strip().upper()

        if seq1 == "" or seq2 == "":
            QMessageBox.critical(self, "Error", "Both sequences are required!")
            return

        sequence_type = "DNA" if self.dna_radio.isChecked() else ("RNA" if self.rna_radio.isChecked() else "Protein")
        alignment_type = "global" if self.global_radio.isChecked() else "local"

        match_penalty = float(self.match_penalty_edit.text())
        mismatch_penalty = float(self.mismatch_penalty_edit.text())
        gap_open_penalty = float(self.gap_open_penalty_edit.text())
        gap_extend_penalty = float(self.gap_extend_penalty_edit.text())

        aligner = PairwiseAligner()
        aligner.mode = alignment_type
        aligner.match_score = match_penalty
        aligner.mismatch_score = mismatch_penalty
        aligner.open_gap_score = gap_open_penalty
        aligner.extend_gap_score = gap_extend_penalty

        if sequence_type == 'DNA':
            valid_chars = set("ATGC")
        elif sequence_type == 'RNA':
            valid_chars = set("AUGC")
        elif sequence_type == 'Protein':
            valid_chars = set("ACDEFGHIKLMNPQRSTVWY")

        if not set(seq1).issubset(valid_chars) or not set(seq2).issubset(valid_chars):
            QMessageBox.critical(self, "Error", f"Invalid characters in sequences! Only {', '.join(valid_chars)} are allowed.")
            return

        alignments = aligner.align(seq1, seq2)

        if not alignments:
            QMessageBox.critical(self, "Error", "No alignment found!")
            return

        max_score = max(alignment.score for alignment in alignments)
        best_alignments = [alignment for alignment in alignments if alignment.score == max_score]

        result_text = ""
        for alignment in best_alignments:
            score_rounded = round(alignment.score, 2)
            result_text += f"<b>Alignment Score:</b> {score_rounded}<br>"
            lines = str(alignment).split('\n')  # Split alignment into lines
            
            for line in lines:
                if not 'target' in line and not 'query' in line:
                    result_text += f"&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<font face='Courier New'>{line}</font><br>"  # Apply tab and Courier New font to each line
                elif 'query' in line:
                    # Split the line at the first occurrence of a digit
                    number_index = next((i for i, c in enumerate(line) if c.isdigit()), None)

                    if number_index is not None:
                        # Add tabs before and after the number
                        line_with_tabs = f"<font face='Courier New'>{line[:number_index]}&nbsp;{line[number_index:]}</font>"
                    else:
                        line_with_tabs = f"<font face='Courier New'>{line}</font>"

                    result_text += f"{line_with_tabs}<br>"  # Apply Courier New font to each line with tabs
                else:
                    result_text += f"<font face='Courier New'>{line}</font><br>"  # Apply Courier New font to each line
            result_text += "<br>"  # Add extra line break after each alignment

        self.result_edit.setHtml(result_text)
        self.result_label.show()
        self.result_edit.show()
        self.save_button.show() 

    def upload_sequence(self, text_edit):
        file_dialog = QFileDialog()
        file_dialog.setNameFilter("Sequence files (*.txt *.fasta)")
        file_dialog.setViewMode(QFileDialog.List)
        file_dialog.setFileMode(QFileDialog.ExistingFile)

        if file_dialog.exec_():
            file_path = file_dialog.selectedFiles()[0]
            with open(file_path, 'r') as file:
                sequence = ''
                found_first_sequence = False
                for line in file:
                    if line.startswith('>'):
                        if not found_first_sequence:
                            found_first_sequence = True
                        else:
                            break
                        continue
                    if found_first_sequence:
                        sequence += line.strip()
                text_edit.setPlainText(sequence)
            
    def reset_fields(self):
        self.seq1_edit.clear()
        self.seq2_edit.clear()
        self.result_edit.clear()
        self.match_penalty_edit.setText("1")
        self.mismatch_penalty_edit.setText("-1")
        self.gap_open_penalty_edit.setText("-0.5")
        self.gap_extend_penalty_edit.setText("-0.1")

        # Reset radio buttons to default
        self.dna_radio.setChecked(True)  # Resets to DNA being selected
        self.global_radio.setChecked(True) 
        
        self.result_label.hide()  # Hide the label initially
        self.result_edit.hide()  # Hide initially
        self.result_edit.hide()  # Hide initially

    def save_results(self):
        file_path, _ = QFileDialog.getSaveFileName(self, "Save File", "", "Text Files (*.txt);;All Files (*)")
        if file_path:
            with open(file_path, 'w') as file:
                file.write(self.result_edit.toPlainText())
                QMessageBox.information(self, "Save Successful", "The alignment results were successfully saved!")

def main():
    app = QApplication(sys.argv)
    ex = SeqAlignApp()
    ex.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
