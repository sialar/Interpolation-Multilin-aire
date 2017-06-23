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

    void mx(int m) { m0 = m; }
    void my(int m) { m1 = m; }
    void mz(int m) { m2 = m; }
    void mt(int m) { m3 = m; }
    void ms(int m) { m4 = m; }

    void updateX(double _x) { x = _x; }
    void updateY(double _y) { y = _y; }
    void updateZ(double _z) { z = _z; }
    void updateT(double _t) { t = _t; }
    void updateS(double _s) { s = _s; }

    void updateN(int _n) { n = _n; }
    void updateD(int _d) { d = _d; updateCoordsField(); }
    void updateMaxIter(int _nIter) { maxIteration = _nIter; }

    void plotErrors();
    void plotResults();
    void randomPoint();
    void startInterpolation();
    void plotInterpolationPoints();

private:
    double x, y, z, t, s;
    int m0, m1, m2, m3, m4;
    int d, n, f, maxIteration;
    string real_f, f_tilde, relative_error, mse_error, exec_time;
    QSpinBox *mxSpinBox, *mySpinBox, *mzSpinBox, *mtSpinBox, *msSpinBox;
    QDoubleSpinBox *xSpinBox, *ySpinBox, *zSpinBox, *tSpinBox, *sSpinBox;
    QTextEdit *approxValue, *exactValue, *rErrorValue, *mseErrorValue, *execTimeValue;


    void updateGUI();
    void updateCoordsField();
    void createControlButtons(int a, int b, int c, int d);
    void createInputFieldForMethod(int a, int b, int c, int d);
    void createOutputFieldForResults(int a, int b, int c, int d);
    void createInputFieldsForParameters(int a, int b, int c, int d);
    void createExclusiveGroupForFunctions(int a, int b, int c, int d);
    void createInputFieldForPointCoordinates(int a, int b, int c, int d);
};

#endif
