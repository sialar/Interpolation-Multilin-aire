#ifndef WINDOW_H
#define WINDOW_H

#include <QWidget>
#include <QApplication>
#include <QWidget>
#include <QDoubleSpinBox>
#include <QPushButton>
#include <QTextEdit>
#include <QSlider>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

class QGroupBox;

class Window : public QWidget
{
    Q_OBJECT

public:
    Window();
    ~Window();
    static std::string projectPath;


public slots:
    void f1() { f = 1; }
    void f2() { f = 2; }
    void f3() { f = 3; }

    void m0() { method = 0; }
    void m1() { method = 1; }
    void m2() { method = 2; }

    void updateX(double _x) { x = _x; }
    void updateY(double _y) { y = _y; }
    void updateZ(double _z) { z = _z; }
    void updateT(double _t) { t = _t; }
    void updateS(double _s) { s = _s; }

    void updateMaxIter(int _nIter) { maxIteration = _nIter; }
    void updateD(int _d) { d = _d; updateCoordsField(); }
    void updateN(int _n) { n = _n; }

    void startInterpolation();
    void plotResults();
    void plotInterpolationPoints();
    void plotErrors();
    void randomPoint();

private:
    int f;
    int method;
    int maxIteration;
    int d, n;
    double x, y, z, t, s;
    QDoubleSpinBox *xSpinBox, *ySpinBox, *zSpinBox, *tSpinBox, *sSpinBox;
    QSpinBox *mxSpinBox, *mySpinBox, *mzSpinBox, *mtSpinBox, *msSpinBox;
    string real_f;
    string f_tilde;
    string relative_error;
    string mse_error;
    string exec_time;

    QTextEdit *approxValue, *exactValue, *rErrorValue, *mseErrorValue, *execTimeValue;

    void updateGUI();
    void updateCoordsField();
    void createControlButtons(int a, int b, int c, int d);
    void createExclusiveGroupForFunctions(int a, int b, int c, int d);
    void createInputFieldForPointCoordinates(int a, int b, int c, int d);
    void createInputFieldForMethod(int a, int b, int c, int d);
    void createInputFieldsForParameters(int a, int b, int c, int d);
    void createOutputFieldForResults(int a, int b, int c, int d);
};

#endif
