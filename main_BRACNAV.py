#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from PyQt5 import QtWidgets,QtGui,QtCore
from PyQt5.QtCore import QThread
##,QObject,pyqtSlot,pyqtSignal
import mainwindow
import bracnav
class BRACNAVmain(QtWidgets.QMainWindow, mainwindow.Ui_BRACNAV):
    def __init__(self):
        super().__init__()
        self.setupUi(self)
        self.show()
        self.input_file_button.clicked.connect(self.input_file_button_onClick)
        self.amplicon_file_button.clicked.connect(self.amplicon_file_button_onClick)
        self.patients_table_button.clicked.connect(self.patients_table_button_onClick)
        self.output_file_button.clicked.connect(self.output_file_button_onClick)
        self.start_button.clicked.connect(self.start_button_onClick)
        self.progressBar.setEnabled(False)
    def input_file_button_onClick(self):
        dialog=QtWidgets.QFileDialog()
        self.input_file_path=dialog.getOpenFileName(self,'Select directory with amplicon coverages files')
        self.input_file_text.setText(self.input_file_path[0])
    def amplicon_file_button_onClick(self):
        dialog=QtWidgets.QFileDialog()
        self.amplicon_file_path=dialog.getOpenFileName(self,'Select file with amplicon coordinates')
        self.amplicon_file_text.setText(self.amplicon_file_path[0])
    def patients_table_button_onClick(self):
        dialog=QtWidgets.QFileDialog()
        self.patients_table_path=dialog.getOpenFileName(self,'Select file with patient table')
        self.patients_table_text.setText(self.patients_table_path[0])
    def output_file_button_onClick(self):
        dialog=QtWidgets.QFileDialog()
        self.output_file_path=dialog.getOpenFileName(self,'Select output file')
        self.output_file_text.setText(self.output_file_path[0])
    def start_button_onClick(self):
        self.thread=QThread()
        self.input_file_button.setEnabled(False)
        self.amplicon_file_button.setEnabled(False)
        self.patients_table_button.setEnabled(False)
        self.output_file_button.setEnabled(False)
        self.start_button.setEnabled(False)
        self.reference_version_box.setEnabled(False)
        self.clust_box.setEnabled(False)
        self.progressBar.setEnabled(True)
        self.bracnavObject=bracnav.BRACNAV(self.input_file_text.text(),
                                           self.amplicon_file_text.text(),
                                           self.patients_table_text.text(),
                                           self.output_file_text.text(),
                                           self.reference_version_box.currentText(),
                                           self.clust_box.isChecked(),
                                           float(self.Hscores_threshold_text.text()),
                                           float(self.Hpvalue_text.text()),
                                           float(self.minimal_score_text.text()),
                                           float(self.maximal_p_value_text.text()),
                                           int(self.minimal_median_coverage_text.text()),
                                           int(self.permutations_number_text.text()),
                                           int(self.exonCoveredWhole_text.text()),
                                           float(self.exonCoveredPart_text.text()),
                                           int(self.exonNonCovered_text.text()),
                                           float(self.delTh1_text.text()),
                                           float(self.delTh2_text.text()),
                                           float(self.duplTh1_text.text()),
                                           float(self.duplTh2_text.text()),
                                           float(self.delta_text.text()))
        self.bracnavObject.moveToThread(self.thread)
        self.thread.started.connect(self.bracnavObject.run)
        self.bracnavObject.finished.connect(self.thread.quit)
        self.bracnavObject.done_percent.connect(self.progressBar.setValue)
        self.bracnavObject.finished.connect(self.bracnavObject.deleteLater)
        self.thread.finished.connect(self.thread.deleteLater)
        self.thread.finished.connect(self.returnStartState)
        print('Starting thread...')
        self.thread.start()
    def returnStartState(self):
        try:
            self.thread.wait()
        except RuntimeError:
            pass
        print('\nReturning start state of buttons and labels...')
        self.input_file_button.setEnabled(True)
        self.amplicon_file_button.setEnabled(True)
        self.patients_table_button.setEnabled(True)
        self.output_file_button.setEnabled(True)
        self.start_button.setEnabled(True)
        self.reference_version_box.setEnabled(True)
        self.clust_box.setEnabled(True)  
        
def main():
    app = QtWidgets.QApplication(sys.argv)
    window = BRACNAVmain()
    app.exec_()

if __name__ == '__main__':
        main()
