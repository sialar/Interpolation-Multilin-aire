#include <QtWidgets>

#include "interface.h"

Window::Window() : QWidget()
{
    setFixedSize(1200, 675);
    createExclusiveGroupForFunctions();
    createExclusiveGroupForMethods();
}

QGroupBox *Window::createExclusiveGroupForFunctions()
{
    QGroupBox *groupBox = new QGroupBox(tr("Choose the function to interpolate at point X=(x,y,z)"), this);
    QRadioButton *radio1 = new QRadioButton(tr("f1(x) = random_polynomial of degree 7"));
    QRadioButton *radio2 = new QRadioButton(tr("f2(x) = sqrt(1-x*x)*exp(-y-z)"));
    QRadioButton *radio3 = new QRadioButton(tr("f3(x) = sqrt(||X||)"));

    connect(radio1, SIGNAL(clicked()), this, SLOT(chooseFunction1()));
    connect(radio2, SIGNAL(clicked()), this, SLOT(chooseFunction2()));
    connect(radio3, SIGNAL(clicked()), this, SLOT(chooseFunction3()));

    radio1->setChecked(true);

    QVBoxLayout *vbox = new QVBoxLayout;
    vbox->addWidget(radio1);
    vbox->addWidget(radio2);
    vbox->addWidget(radio3);
    vbox->addStretch(1);
    groupBox->setLayout(vbox);

    groupBox->setFixedSize(320,100);
    groupBox->move(50,50);
    return groupBox;
}

QGroupBox *Window::createExclusiveGroupForMethods()
{
    QGroupBox *groupBox = new QGroupBox(tr("Choose the method of interpolation"), this);
    QRadioButton *radio1 = new QRadioButton(tr("M0 = Using lagrange polynomial functions and leja points"));
    QRadioButton *radio2 = new QRadioButton(tr("M1 = Using piecewise linear functions and middle points"));
    QRadioButton *radio3 = new QRadioButton(tr("M2 = Using piecewise quadratic functions and middle points"));

    connect(radio1, SIGNAL(clicked()), this, SLOT(chooseMethod0()));
    connect(radio2, SIGNAL(clicked()), this, SLOT(chooseMethod1()));
    connect(radio3, SIGNAL(clicked()), this, SLOT(chooseMethod2()));

    radio1->setChecked(true);

    QVBoxLayout *vbox = new QVBoxLayout;
    vbox->addWidget(radio1);
    vbox->addWidget(radio2);
    vbox->addWidget(radio3);
    vbox->addStretch(1);
    groupBox->setLayout(vbox);

    groupBox->setFixedSize(320,100);
    groupBox->move(420,50);
    return groupBox;
}

void Window::chooseFunction1()
{
    std::ofstream out("../gui/txt");
    std::cout << "1" << std::endl;
    out << "1" << std::endl;
}
void Window::chooseFunction2()
{
    std::ofstream out("../gui/txt");
    std::cout << "2" << std::endl;
    out << "2" << std::endl;
}
void Window::chooseFunction3()
{
    std::ofstream out("../gui/txt");
    std::cout << "3" << std::endl;
    out << "3" << std::endl;
}
void Window::chooseMethod0()
{
    std::ofstream out("../gui/txt");
    std::cout << "0" << std::endl;
    out << "0" << std::endl;
}
void Window::chooseMethod1()
{
    std::ofstream out("../gui/txt");
    std::cout << "1" << std::endl;
    out << "1" << std::endl;
}
void Window::chooseMethod2()
{
    std::ofstream out("../gui/txt");
    std::cout << "2" << std::endl;
    out << "2" << std::endl;
}

Window::~Window()
{

}
