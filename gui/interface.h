#ifndef WINDOW_H
#define WINDOW_H

#include <QWidget>
#include <QApplication>
#include <QWidget>
#include <QPushButton>
#include <QSlider>
#include <iostream>
#include <fstream>


class QGroupBox;

class Window : public QWidget
{
    Q_OBJECT

public:
    Window();
    ~Window();

public slots:
    void chooseFunction1();
    void chooseFunction2();
    void chooseFunction3();
    void chooseMethod0();
    void chooseMethod1();
    void chooseMethod2();

private:
    QGroupBox *createExclusiveGroupForFunctions();
    QGroupBox *createExclusiveGroupForMethods();
};

#endif
