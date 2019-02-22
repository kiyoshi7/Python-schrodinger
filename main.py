# ------------------------------------------------------
# ---------------------- main.py -----------------------
# ------------------------------------------------------
from PyQt5.QtWidgets import*
from PyQt5.uic import loadUi

from matplotlib.backends.backend_qt5agg import (NavigationToolbar2QT as NavigationToolbar)

import numpy as np
import random

## User Imports ##
from Schroedinger_solveIVP import TimeIndependent
     
class MatplotlibWidget(QMainWindow):
    PotentialFunction = 0
    
    def __init__(self):
        self.m = TimeIndependent()
        QMainWindow.__init__(self)

        loadUi("UI.ui",self)

        self.setWindowTitle("PyQt5 & Matplotlib Example GUI")

        #self.pushButton_generate_random_signal.clicked.connect(self.PlotPotentialFunction)
        ## StepFunction Values
        self.BoxVmax.setPlainText(str(self.m.Vmax))
        self.BoxVmin.setPlainText(str(self.m.Vmin))
        self.BoxVo.setPlainText(str(self.m.Vb))
        self.BoxLenA.setPlainText(str(self.m.SA))
        self.BoxLenB.setPlainText(str(self.m.SB))
        self.BoxSRange.setPlainText(str(self.m.Range))
        self.PlotStepFunction.clicked.connect(lambda: self.PlotPotentialFunction(1))

        ## Quadratic Values
        self.BoxQB4.setPlainText(str(self.m.QB4))
        self.BoxQB3.setPlainText(str(self.m.QB3))
        self.BoxQB2.setPlainText(str(self.m.QB2))
        self.BoxQB1.setPlainText(str(self.m.QB1))
        self.BoxQB0.setPlainText(str(self.m.QB0))
        self.BoxQXRange.setPlainText(str(self.m.Range))
        self.BoxEMax.setPlainText(str(self.m.Vmax))
        self.BoxEMin.setPlainText(str(self.m.Vmin))
        self.PlotQuadraticFunction.clicked.connect(lambda: self.PlotPotentialFunction(2))

        ## Results
        self.ButtonCalculate.clicked.connect(self.PlotSimulation)
        
        
        self.addToolBar(NavigationToolbar(self.MplWidget.canvas, self))
        self.PlotPotentialFunction(0)

    def PlotPotentialFunction(self, Option):
        self.MplWidget.canvas.axes.clear()
        self.PotentialFunction = Option
        ## update values for solver
        if Option == 1:
            
            #self.m.option = 1 # updates the option for SE solver
            self.m.Vmax = float(self.BoxVmax.toPlainText()) ## updates all the constants
            self.m.Vmin = float(self.BoxVmin.toPlainText())
            self.m.Vo   = float(self.BoxVo.toPlainText())
            self.m.SA   = float(self.BoxLenA.toPlainText())
            self.m.SB   = float(self.BoxLenB.toPlainText())
            self.m.Range = float(self.BoxSRange.toPlainText())
            self.MplWidget.canvas.axes.plot(self.m.PotentialStep()[0], self.m.PotentialStep()[1])
            self.MplWidget.canvas.axes.legend(('Step'),loc='upper right')
            
        if Option == 2:
            #self.m.option = 2 # updates the option for SE solver
            self.m.QB4 = float(self.BoxQB4.toPlainText())
            self.m.QB3 = float(self.BoxQB3.toPlainText())
            self.m.QB2 = float(self.BoxQB2.toPlainText())
            self.m.QB1 = float(self.BoxQB1.toPlainText())
            self.m.QB0 = float(self.BoxQB0.toPlainText())
            self.m.Vmax = float(self.BoxEMax.toPlainText())
            self.m.Vmin = float(self.BoxEMin.toPlainText())
            self.m.Vb = float(self.BoxEMin.toPlainText())
            self.m.Range = float(self.BoxQXRange.toPlainText())
            
            self.MplWidget.canvas.axes.plot(self.m.QuadraticFunction()[0], self.m.QuadraticFunction()[1])
            self.MplWidget.canvas.axes.legend(('Quadratic'),loc='upper right')
        self.MplWidget.canvas.draw()
        
    def PlotSimulation(self):
        self.MplWidget.canvas.axes.clear()
        self.PlotPotentialFunction(self.PotentialFunction)

        if self.PotentialFunction == 0:
            self.tabWidget.setCurrentIndex(0)
        else:
            self.tabWidget.setCurrentIndex(3)
            results = self.m.calculate()
            x = results[0]
            y = results[1]
            for i in range(len(x)):
                self.MplWidget.canvas.axes.plot(x[i],y[i])
            self.MplWidget.canvas.draw()
            EnergiesFound = QTextEdit()
            self.ContainerSimulated.setWidget(EnergiesFound)
            for i in range(len(x)):
                EnergiesFound.append(str(y[i][0]))

        
        

app = QApplication([])
window = MatplotlibWidget()
window.show()
app.exec_()
