# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'mainwindow.ui'
#
# Created by: PyQt5 UI code generator 5.15.2
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_BRACNAV(object):
    def setupUi(self, BRACNAV):
        BRACNAV.setObjectName("BRACNAV")
        BRACNAV.resize(800, 404)
        BRACNAV.setMaximumSize(QtCore.QSize(800, 404))
        self.centralWidget = QtWidgets.QWidget(BRACNAV)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.centralWidget.sizePolicy().hasHeightForWidth())
        self.centralWidget.setSizePolicy(sizePolicy)
        self.centralWidget.setObjectName("centralWidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.centralWidget)
        self.verticalLayout.setContentsMargins(11, 11, 11, 11)
        self.verticalLayout.setSpacing(6)
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setSpacing(6)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.input_file_text = QtWidgets.QLineEdit(self.centralWidget)
        self.input_file_text.setEnabled(True)
        self.input_file_text.setMinimumSize(QtCore.QSize(200, 25))
        self.input_file_text.setMaximumSize(QtCore.QSize(300, 20))
        self.input_file_text.setObjectName("input_file_text")
        self.horizontalLayout.addWidget(self.input_file_text)
        self.label = QtWidgets.QLabel(self.centralWidget)
        self.label.setMaximumSize(QtCore.QSize(300, 30))
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.label.setFont(font)
        self.label.setObjectName("label")
        self.horizontalLayout.addWidget(self.label)
        self.input_file_button = QtWidgets.QPushButton(self.centralWidget)
        self.input_file_button.setMaximumSize(QtCore.QSize(120, 30))
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.input_file_button.setFont(font)
        self.input_file_button.setObjectName("input_file_button")
        self.horizontalLayout.addWidget(self.input_file_button)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setSpacing(6)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.amplicon_file_text = QtWidgets.QLineEdit(self.centralWidget)
        self.amplicon_file_text.setEnabled(True)
        self.amplicon_file_text.setMinimumSize(QtCore.QSize(200, 25))
        self.amplicon_file_text.setMaximumSize(QtCore.QSize(300, 20))
        self.amplicon_file_text.setObjectName("amplicon_file_text")
        self.horizontalLayout_2.addWidget(self.amplicon_file_text)
        self.label_2 = QtWidgets.QLabel(self.centralWidget)
        self.label_2.setMaximumSize(QtCore.QSize(300, 30))
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.label_2.setFont(font)
        self.label_2.setScaledContents(False)
        self.label_2.setObjectName("label_2")
        self.horizontalLayout_2.addWidget(self.label_2)
        self.amplicon_file_button = QtWidgets.QPushButton(self.centralWidget)
        self.amplicon_file_button.setMaximumSize(QtCore.QSize(120, 16777215))
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.amplicon_file_button.setFont(font)
        self.amplicon_file_button.setObjectName("amplicon_file_button")
        self.horizontalLayout_2.addWidget(self.amplicon_file_button)
        self.verticalLayout.addLayout(self.horizontalLayout_2)
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setContentsMargins(-1, 0, -1, -1)
        self.horizontalLayout_3.setSpacing(6)
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.patients_table_text = QtWidgets.QLineEdit(self.centralWidget)
        self.patients_table_text.setEnabled(True)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.patients_table_text.sizePolicy().hasHeightForWidth())
        self.patients_table_text.setSizePolicy(sizePolicy)
        self.patients_table_text.setMinimumSize(QtCore.QSize(0, 25))
        self.patients_table_text.setMaximumSize(QtCore.QSize(300, 25))
        self.patients_table_text.setObjectName("patients_table_text")
        self.horizontalLayout_3.addWidget(self.patients_table_text)
        self.label_3 = QtWidgets.QLabel(self.centralWidget)
        self.label_3.setMaximumSize(QtCore.QSize(300, 30))
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.label_3.setFont(font)
        self.label_3.setObjectName("label_3")
        self.horizontalLayout_3.addWidget(self.label_3)
        self.patients_table_button = QtWidgets.QPushButton(self.centralWidget)
        self.patients_table_button.setMaximumSize(QtCore.QSize(120, 30))
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.patients_table_button.setFont(font)
        self.patients_table_button.setObjectName("patients_table_button")
        self.horizontalLayout_3.addWidget(self.patients_table_button)
        self.verticalLayout.addLayout(self.horizontalLayout_3)
        self.horizontalLayout_4 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_4.setContentsMargins(-1, -1, -1, 0)
        self.horizontalLayout_4.setSpacing(6)
        self.horizontalLayout_4.setObjectName("horizontalLayout_4")
        self.output_file_text = QtWidgets.QLineEdit(self.centralWidget)
        self.output_file_text.setEnabled(True)
        self.output_file_text.setMaximumSize(QtCore.QSize(300, 25))
        self.output_file_text.setBaseSize(QtCore.QSize(220, 0))
        self.output_file_text.setObjectName("output_file_text")
        self.horizontalLayout_4.addWidget(self.output_file_text)
        self.label_4 = QtWidgets.QLabel(self.centralWidget)
        self.label_4.setMaximumSize(QtCore.QSize(300, 30))
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.label_4.setFont(font)
        self.label_4.setObjectName("label_4")
        self.horizontalLayout_4.addWidget(self.label_4)
        self.output_file_button = QtWidgets.QPushButton(self.centralWidget)
        self.output_file_button.setMaximumSize(QtCore.QSize(120, 30))
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.output_file_button.setFont(font)
        self.output_file_button.setObjectName("output_file_button")
        self.horizontalLayout_4.addWidget(self.output_file_button)
        self.verticalLayout.addLayout(self.horizontalLayout_4)
        self.horizontalLayout_5 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_5.setContentsMargins(-1, -1, -1, 0)
        self.horizontalLayout_5.setSpacing(6)
        self.horizontalLayout_5.setObjectName("horizontalLayout_5")
        self.reference_version_box = QtWidgets.QComboBox(self.centralWidget)
        self.reference_version_box.setMaximumSize(QtCore.QSize(100, 100))
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.reference_version_box.setFont(font)
        self.reference_version_box.setObjectName("reference_version_box")
        self.reference_version_box.addItem("")
        self.reference_version_box.addItem("")
        self.horizontalLayout_5.addWidget(self.reference_version_box)
        self.label_5 = QtWidgets.QLabel(self.centralWidget)
        self.label_5.setMaximumSize(QtCore.QSize(200, 30))
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.label_5.setFont(font)
        self.label_5.setObjectName("label_5")
        self.horizontalLayout_5.addWidget(self.label_5)
        self.clust_box = QtWidgets.QCheckBox(self.centralWidget)
        self.clust_box.setMaximumSize(QtCore.QSize(200, 30))
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.clust_box.setFont(font)
        self.clust_box.setObjectName("clust_box")
        self.horizontalLayout_5.addWidget(self.clust_box)
        self.verticalLayout.addLayout(self.horizontalLayout_5)
        self.horizontalLayout_6 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_6.setContentsMargins(-1, -1, -1, 0)
        self.horizontalLayout_6.setSpacing(0)
        self.horizontalLayout_6.setObjectName("horizontalLayout_6")
        self.Hscores_threshold_text = QtWidgets.QLineEdit(self.centralWidget)
        self.Hscores_threshold_text.setMaximumSize(QtCore.QSize(65, 25))
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.Hscores_threshold_text.setFont(font)
        self.Hscores_threshold_text.setAlignment(QtCore.Qt.AlignCenter)
        self.Hscores_threshold_text.setObjectName("Hscores_threshold_text")
        self.horizontalLayout_6.addWidget(self.Hscores_threshold_text)
        self.label_7 = QtWidgets.QLabel(self.centralWidget)
        self.label_7.setMaximumSize(QtCore.QSize(155, 16777215))
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.label_7.setFont(font)
        self.label_7.setObjectName("label_7")
        self.horizontalLayout_6.addWidget(self.label_7)
        self.Hpvalue_text = QtWidgets.QLineEdit(self.centralWidget)
        self.Hpvalue_text.setMaximumSize(QtCore.QSize(65, 25))
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.Hpvalue_text.setFont(font)
        self.Hpvalue_text.setAlignment(QtCore.Qt.AlignCenter)
        self.Hpvalue_text.setObjectName("Hpvalue_text")
        self.horizontalLayout_6.addWidget(self.Hpvalue_text)
        self.label_6 = QtWidgets.QLabel(self.centralWidget)
        self.label_6.setMaximumSize(QtCore.QSize(160, 16777215))
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.label_6.setFont(font)
        self.label_6.setObjectName("label_6")
        self.horizontalLayout_6.addWidget(self.label_6)
        self.permutations_number_text = QtWidgets.QLineEdit(self.centralWidget)
        self.permutations_number_text.setEnabled(True)
        self.permutations_number_text.setMaximumSize(QtCore.QSize(65, 25))
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.permutations_number_text.setFont(font)
        self.permutations_number_text.setAlignment(QtCore.Qt.AlignCenter)
        self.permutations_number_text.setObjectName("permutations_number_text")
        self.horizontalLayout_6.addWidget(self.permutations_number_text)
        self.label_8 = QtWidgets.QLabel(self.centralWidget)
        self.label_8.setMaximumSize(QtCore.QSize(200, 100))
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.label_8.setFont(font)
        self.label_8.setObjectName("label_8")
        self.horizontalLayout_6.addWidget(self.label_8)
        self.verticalLayout.addLayout(self.horizontalLayout_6)
        self.horizontalLayout_7 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_7.setContentsMargins(-1, -1, -1, 0)
        self.horizontalLayout_7.setSpacing(0)
        self.horizontalLayout_7.setObjectName("horizontalLayout_7")
        self.minimal_score_text = QtWidgets.QLineEdit(self.centralWidget)
        self.minimal_score_text.setMaximumSize(QtCore.QSize(65, 25))
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.minimal_score_text.setFont(font)
        self.minimal_score_text.setAlignment(QtCore.Qt.AlignCenter)
        self.minimal_score_text.setObjectName("minimal_score_text")
        self.horizontalLayout_7.addWidget(self.minimal_score_text)
        self.label_9 = QtWidgets.QLabel(self.centralWidget)
        self.label_9.setMaximumSize(QtCore.QSize(155, 16777215))
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.label_9.setFont(font)
        self.label_9.setObjectName("label_9")
        self.horizontalLayout_7.addWidget(self.label_9)
        self.maximal_p_value_text = QtWidgets.QLineEdit(self.centralWidget)
        self.maximal_p_value_text.setMaximumSize(QtCore.QSize(65, 25))
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.maximal_p_value_text.setFont(font)
        self.maximal_p_value_text.setAlignment(QtCore.Qt.AlignCenter)
        self.maximal_p_value_text.setObjectName("maximal_p_value_text")
        self.horizontalLayout_7.addWidget(self.maximal_p_value_text)
        self.label_10 = QtWidgets.QLabel(self.centralWidget)
        self.label_10.setMaximumSize(QtCore.QSize(160, 16777215))
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.label_10.setFont(font)
        self.label_10.setObjectName("label_10")
        self.horizontalLayout_7.addWidget(self.label_10)
        self.minimal_median_coverage_text = QtWidgets.QLineEdit(self.centralWidget)
        self.minimal_median_coverage_text.setMaximumSize(QtCore.QSize(65, 25))
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.minimal_median_coverage_text.setFont(font)
        self.minimal_median_coverage_text.setAlignment(QtCore.Qt.AlignCenter)
        self.minimal_median_coverage_text.setObjectName("minimal_median_coverage_text")
        self.horizontalLayout_7.addWidget(self.minimal_median_coverage_text)
        self.label_11 = QtWidgets.QLabel(self.centralWidget)
        self.label_11.setMaximumSize(QtCore.QSize(200, 30))
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.label_11.setFont(font)
        self.label_11.setObjectName("label_11")
        self.horizontalLayout_7.addWidget(self.label_11)
        self.verticalLayout.addLayout(self.horizontalLayout_7)
        self.horizontalLayout_10 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_10.setContentsMargins(-1, -1, -1, 0)
        self.horizontalLayout_10.setSpacing(0)
        self.horizontalLayout_10.setObjectName("horizontalLayout_10")
        self.exonCoveredWhole_text = QtWidgets.QLineEdit(self.centralWidget)
        self.exonCoveredWhole_text.setMaximumSize(QtCore.QSize(65, 25))
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.exonCoveredWhole_text.setFont(font)
        self.exonCoveredWhole_text.setAlignment(QtCore.Qt.AlignCenter)
        self.exonCoveredWhole_text.setObjectName("exonCoveredWhole_text")
        self.horizontalLayout_10.addWidget(self.exonCoveredWhole_text)
        self.label_14 = QtWidgets.QLabel(self.centralWidget)
        self.label_14.setMaximumSize(QtCore.QSize(155, 16777215))
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.label_14.setFont(font)
        self.label_14.setObjectName("label_14")
        self.horizontalLayout_10.addWidget(self.label_14)
        self.exonCoveredPart_text = QtWidgets.QLineEdit(self.centralWidget)
        self.exonCoveredPart_text.setMaximumSize(QtCore.QSize(65, 25))
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.exonCoveredPart_text.setFont(font)
        self.exonCoveredPart_text.setAlignment(QtCore.Qt.AlignCenter)
        self.exonCoveredPart_text.setObjectName("exonCoveredPart_text")
        self.horizontalLayout_10.addWidget(self.exonCoveredPart_text)
        self.label_13 = QtWidgets.QLabel(self.centralWidget)
        self.label_13.setMaximumSize(QtCore.QSize(160, 16777215))
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.label_13.setFont(font)
        self.label_13.setObjectName("label_13")
        self.horizontalLayout_10.addWidget(self.label_13)
        self.exonNonCovered_text = QtWidgets.QLineEdit(self.centralWidget)
        self.exonNonCovered_text.setMaximumSize(QtCore.QSize(65, 25))
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.exonNonCovered_text.setFont(font)
        self.exonNonCovered_text.setAlignment(QtCore.Qt.AlignCenter)
        self.exonNonCovered_text.setObjectName("exonNonCovered_text")
        self.horizontalLayout_10.addWidget(self.exonNonCovered_text)
        self.label_12 = QtWidgets.QLabel(self.centralWidget)
        self.label_12.setMaximumSize(QtCore.QSize(200, 16777215))
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.label_12.setFont(font)
        self.label_12.setObjectName("label_12")
        self.horizontalLayout_10.addWidget(self.label_12)
        self.verticalLayout.addLayout(self.horizontalLayout_10)
        self.horizontalLayout_11 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_11.setContentsMargins(-1, -1, -1, 0)
        self.horizontalLayout_11.setSpacing(0)
        self.horizontalLayout_11.setObjectName("horizontalLayout_11")
        self.delTh1_text = QtWidgets.QLineEdit(self.centralWidget)
        self.delTh1_text.setMaximumSize(QtCore.QSize(65, 25))
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.delTh1_text.setFont(font)
        self.delTh1_text.setAlignment(QtCore.Qt.AlignCenter)
        self.delTh1_text.setObjectName("delTh1_text")
        self.horizontalLayout_11.addWidget(self.delTh1_text)
        self.label_15 = QtWidgets.QLabel(self.centralWidget)
        self.label_15.setMaximumSize(QtCore.QSize(65, 30))
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.label_15.setFont(font)
        self.label_15.setObjectName("label_15")
        self.horizontalLayout_11.addWidget(self.label_15)
        self.delTh2_text = QtWidgets.QLineEdit(self.centralWidget)
        self.delTh2_text.setMaximumSize(QtCore.QSize(65, 25))
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.delTh2_text.setFont(font)
        self.delTh2_text.setAlignment(QtCore.Qt.AlignCenter)
        self.delTh2_text.setObjectName("delTh2_text")
        self.horizontalLayout_11.addWidget(self.delTh2_text)
        self.label_16 = QtWidgets.QLabel(self.centralWidget)
        self.label_16.setMaximumSize(QtCore.QSize(65, 16777215))
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.label_16.setFont(font)
        self.label_16.setObjectName("label_16")
        self.horizontalLayout_11.addWidget(self.label_16)
        self.duplTh1_text = QtWidgets.QLineEdit(self.centralWidget)
        self.duplTh1_text.setMaximumSize(QtCore.QSize(65, 25))
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.duplTh1_text.setFont(font)
        self.duplTh1_text.setAlignment(QtCore.Qt.AlignCenter)
        self.duplTh1_text.setObjectName("duplTh1_text")
        self.horizontalLayout_11.addWidget(self.duplTh1_text)
        self.label_17 = QtWidgets.QLabel(self.centralWidget)
        self.label_17.setMaximumSize(QtCore.QSize(65, 16777215))
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.label_17.setFont(font)
        self.label_17.setObjectName("label_17")
        self.horizontalLayout_11.addWidget(self.label_17)
        self.duplTh2_text = QtWidgets.QLineEdit(self.centralWidget)
        self.duplTh2_text.setMaximumSize(QtCore.QSize(65, 25))
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.duplTh2_text.setFont(font)
        self.duplTh2_text.setAlignment(QtCore.Qt.AlignCenter)
        self.duplTh2_text.setObjectName("duplTh2_text")
        self.horizontalLayout_11.addWidget(self.duplTh2_text)
        self.label_18 = QtWidgets.QLabel(self.centralWidget)
        self.label_18.setMaximumSize(QtCore.QSize(65, 16777215))
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.label_18.setFont(font)
        self.label_18.setObjectName("label_18")
        self.horizontalLayout_11.addWidget(self.label_18)
        self.delta_text = QtWidgets.QLineEdit(self.centralWidget)
        self.delta_text.setMaximumSize(QtCore.QSize(65, 25))
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.delta_text.setFont(font)
        self.delta_text.setAlignment(QtCore.Qt.AlignCenter)
        self.delta_text.setObjectName("delta_text")
        self.horizontalLayout_11.addWidget(self.delta_text)
        self.delta = QtWidgets.QLabel(self.centralWidget)
        self.delta.setMaximumSize(QtCore.QSize(65, 16777215))
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.delta.setFont(font)
        self.delta.setObjectName("delta")
        self.horizontalLayout_11.addWidget(self.delta)
        self.verticalLayout.addLayout(self.horizontalLayout_11)
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem)
        self.horizontalLayout_9 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_9.setContentsMargins(12, -1, 12, 0)
        self.horizontalLayout_9.setSpacing(6)
        self.horizontalLayout_9.setObjectName("horizontalLayout_9")
        self.progressBar = QtWidgets.QProgressBar(self.centralWidget)
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.progressBar.setFont(font)
        self.progressBar.setProperty("value", 0)
        self.progressBar.setObjectName("progressBar")
        self.horizontalLayout_9.addWidget(self.progressBar)
        self.verticalLayout.addLayout(self.horizontalLayout_9)
        spacerItem1 = QtWidgets.QSpacerItem(20, 30, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem1)
        self.horizontalLayout_8 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_8.setContentsMargins(-1, -1, -1, 0)
        self.horizontalLayout_8.setSpacing(6)
        self.horizontalLayout_8.setObjectName("horizontalLayout_8")
        spacerItem2 = QtWidgets.QSpacerItem(300, 0, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_8.addItem(spacerItem2)
        self.start_button = QtWidgets.QPushButton(self.centralWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.start_button.sizePolicy().hasHeightForWidth())
        self.start_button.setSizePolicy(sizePolicy)
        self.start_button.setMaximumSize(QtCore.QSize(250, 30))
        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        self.start_button.setFont(font)
        self.start_button.setObjectName("start_button")
        self.horizontalLayout_8.addWidget(self.start_button)
        spacerItem3 = QtWidgets.QSpacerItem(300, 30, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_8.addItem(spacerItem3)
        self.verticalLayout.addLayout(self.horizontalLayout_8)
        BRACNAV.setCentralWidget(self.centralWidget)
        self.statusBar = QtWidgets.QStatusBar(BRACNAV)
        self.statusBar.setObjectName("statusBar")
        BRACNAV.setStatusBar(self.statusBar)

        self.retranslateUi(BRACNAV)
        QtCore.QMetaObject.connectSlotsByName(BRACNAV)

    def retranslateUi(self, BRACNAV):
        _translate = QtCore.QCoreApplication.translate
        BRACNAV.setWindowTitle(_translate("BRACNAV", "BRACNAV"))
        self.input_file_text.setProperty("html", _translate("BRACNAV", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Ubuntu\'; font-size:11pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.label.setToolTip(_translate("BRACNAV", "<html><head/><body><p>File with amplicon coverages with following columns: Patient#, Patient_ID, Barcodes, Median_Coverage, Number_&lt;30, amplicon#1... (Patient_ID, Barcodes, Median_Coverage, Number_&lt;30 are optional columns)</p></body></html>"))
        self.label.setWhatsThis(_translate("BRACNAV", "<html><head/><body><p>File with amplicon coverages</p><p><br/></p></body></html>"))
        self.label.setText(_translate("BRACNAV", "Input file  (.tsv ) *"))
        self.input_file_button.setText(_translate("BRACNAV", "Choose file"))
        self.amplicon_file_text.setProperty("html", _translate("BRACNAV", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Ubuntu\'; font-size:11pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.label_2.setToolTip(_translate("BRACNAV", "<html><head/><body><p>TSV-file with coordinates for each amplicon and numbers of multiplexes. First column is a number of amplicon, second - number of chromosome, next two - coordinates, and last - number of multiplex.</p></body></html>"))
        self.label_2.setWhatsThis(_translate("BRACNAV", "<html><head/><body><p>TSV-file with coordinates for each amplicon and numbers of multiplexes</p></body></html>"))
        self.label_2.setText(_translate("BRACNAV", "Amplicon coordinates (.tsv/ .bed) *"))
        self.amplicon_file_button.setText(_translate("BRACNAV", "Choose file"))
        self.label_3.setToolTip(_translate("BRACNAV", "<html><head/><body><p>Optional TSV-file with patient IDs and names. Column names are required (Patient_num, Patient_ID, Barcode_1, Barcode_2)</p></body></html>"))
        self.label_3.setWhatsThis(_translate("BRACNAV", "<html><head/><body><p>TSV-file with patient IDs and names. Column names are required</p></body></html>"))
        self.label_3.setText(_translate("BRACNAV", "Patients table ( .tsv/ .xls)"))
        self.patients_table_button.setText(_translate("BRACNAV", "Choose file"))
        self.output_file_text.setProperty("html", _translate("BRACNAV", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Ubuntu\'; font-size:11pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.label_4.setText(_translate("BRACNAV", "Output file (.xls)"))
        self.output_file_button.setText(_translate("BRACNAV", "Choose file"))
        self.reference_version_box.setItemText(0, _translate("BRACNAV", "hg19"))
        self.reference_version_box.setItemText(1, _translate("BRACNAV", "hg38"))
        self.label_5.setToolTip(_translate("BRACNAV", "<html><head/><body><p>Version of the reference genome</p></body></html>"))
        self.label_5.setWhatsThis(_translate("BRACNAV", "<html><head/><body><p>Version of the reference genome</p></body></html>"))
        self.label_5.setText(_translate("BRACNAV", "Reference version"))
        self.clust_box.setToolTip(_translate("BRACNAV", "<html><head/><body><p>Use this parameter if you do not want to clusterize samples during normalization</p></body></html>"))
        self.clust_box.setWhatsThis(_translate("BRACNAV", "<html><head/><body><p>Use this parameter if you do not want to clusterize samples during normalization</p></body></html>"))
        self.clust_box.setText(_translate("BRACNAV", "Not clust"))
        self.Hscores_threshold_text.setText(_translate("BRACNAV", "9.9"))
        self.label_7.setToolTip(_translate("BRACNAV", "<html><head/><body><p>Minimal hard score for large rearrangement detection</p></body></html>"))
        self.label_7.setWhatsThis(_translate("BRACNAV", "<html><head/><body><p>Minimal hard score for large rearrangement detection</p></body></html>"))
        self.label_7.setText(_translate("BRACNAV", "Hard scores threshold"))
        self.Hpvalue_text.setText(_translate("BRACNAV", "0.001"))
        self.label_6.setToolTip(_translate("BRACNAV", "<html><head/><body><p>Maximal p-value for hard filtering mutations</p></body></html>"))
        self.label_6.setWhatsThis(_translate("BRACNAV", "<html><head/><body><p>Maximal p-value for hard filtering mutations</p></body></html>"))
        self.label_6.setText(_translate("BRACNAV", "Hard p-value"))
        self.permutations_number_text.setText(_translate("BRACNAV", "1000"))
        self.label_8.setToolTip(_translate("BRACNAV", "<html><head/><body><p>Number of permutations for calculating p-value</p></body></html>"))
        self.label_8.setWhatsThis(_translate("BRACNAV", "<html><head/><body><p>Number of permutations for calculating p-value</p></body></html>"))
        self.label_8.setText(_translate("BRACNAV", "Permutatious number"))
        self.minimal_score_text.setText(_translate("BRACNAV", "2"))
        self.label_9.setToolTip(_translate("BRACNAV", "<html><head/><body><p>Minimal score for large rearrangement detection</p></body></html>"))
        self.label_9.setWhatsThis(_translate("BRACNAV", "<html><head/><body><p>Minimal score for large rearrangement detection</p></body></html>"))
        self.label_9.setText(_translate("BRACNAV", "Minimal score"))
        self.maximal_p_value_text.setText(_translate("BRACNAV", "0.05"))
        self.label_10.setToolTip(_translate("BRACNAV", "<html><head/><body><p>Maximal p-value for filtering mutations</p></body></html>"))
        self.label_10.setWhatsThis(_translate("BRACNAV", "<html><head/><body><p>Maximal p-value for filtering mutations</p></body></html>"))
        self.label_10.setText(_translate("BRACNAV", "Maximal p-value"))
        self.minimal_median_coverage_text.setText(_translate("BRACNAV", "100"))
        self.label_11.setToolTip(_translate("BRACNAV", "<html><head/><body><p>Minimal median coverage for patient</p></body></html>"))
        self.label_11.setWhatsThis(_translate("BRACNAV", "<html><head/><body><p>Minimal median coverage for patient</p></body></html>"))
        self.label_11.setText(_translate("BRACNAV", "Minimal median coverage"))
        self.exonCoveredWhole_text.setText(_translate("BRACNAV", "1"))
        self.label_14.setToolTip(_translate("BRACNAV", "<html><head/><body><p>Value for calculating score when whole exon is covered </p></body></html>"))
        self.label_14.setText(_translate("BRACNAV", "Whole exon covered"))
        self.exonCoveredPart_text.setText(_translate("BRACNAV", "0.4"))
        self.label_13.setToolTip(_translate("BRACNAV", "<html><head/><body><p>Value for calculating score when part of exon is covered</p></body></html>"))
        self.label_13.setText(_translate("BRACNAV", "Exon covered partially"))
        self.exonNonCovered_text.setText(_translate("BRACNAV", "0"))
        self.label_12.setToolTip(_translate("BRACNAV", "<html><head/><body><p>Value for calculating score when exon is not covered</p></body></html>"))
        self.label_12.setText(_translate("BRACNAV", "Non-covered exon"))
        self.delTh1_text.setText(_translate("BRACNAV", "1.3"))
        self.label_15.setToolTip(_translate("BRACNAV", "<html><head/><body><p>Normalized value of coverage for considering an amplicon as likely deleted</p></body></html>"))
        self.label_15.setText(_translate("BRACNAV", "del1"))
        self.delTh2_text.setText(_translate("BRACNAV", "1.8"))
        self.label_16.setToolTip(_translate("BRACNAV", "<html><head/><body><p>Normalized value of coverage for considering an amplicon as probably deleted </p></body></html>"))
        self.label_16.setText(_translate("BRACNAV", "del2"))
        self.duplTh1_text.setText(_translate("BRACNAV", "2.7"))
        self.label_17.setToolTip(_translate("BRACNAV", "<html><head/><body><p>Normalized value of coverage for considering an amplicon as likely duplicated</p></body></html>"))
        self.label_17.setText(_translate("BRACNAV", "dupl1"))
        self.duplTh2_text.setText(_translate("BRACNAV", "2.4"))
        self.label_18.setToolTip(_translate("BRACNAV", "<html><head/><body><p>Normalized value of coverage for considering an amplicon as probably duplicated</p></body></html>"))
        self.label_18.setText(_translate("BRACNAV", "dupl2"))
        self.delta_text.setText(_translate("BRACNAV", "0.05"))
        self.delta.setToolTip(_translate("BRACNAV", "<html><head/><body><p>Minimal relative difference between two values to be considered as significant</p></body></html>"))
        self.delta.setText(_translate("BRACNAV", "delta"))
        self.start_button.setText(_translate("BRACNAV", "Start"))
