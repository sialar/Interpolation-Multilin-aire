#include <QtWidgets>

#include "interface.h"

Window::Window() : QWidget()
{
    setFixedSize(890, 650);
    createExclusiveGroupForFunctions();
    createExclusiveGroupForMethods();
    createInputFieldForPointCoordinates();
    createInputFieldsForParameters();
    createOutputFieldForResults();
}

void Window::createOutputFieldForResults()
{
    QGroupBox *valuesGroupBox = new QGroupBox(tr("Results"), this);
    QGroupBox *errorGroupBox = new QGroupBox(tr("Interpolation Error"), this);

    QLabel *exactValueLabel = new QLabel(tr("Exact value of f(X)"));
    QTextEdit *exactValue = new QTextEdit;
    exactValue->setReadOnly(true);

    QLabel *approxValueLabel = new QLabel(tr("Approx value of f(X)"));
    QTextEdit *approxValue = new QTextEdit;
    approxValue->setReadOnly(true);

    QLabel *rErrorValueLabel = new QLabel(tr("Relative interpolation error"));
    QTextEdit *rErrorValue = new QTextEdit;
    rErrorValue->setReadOnly(true);

    QLabel *mseValueLabel = new QLabel(tr("MSE interpolation error"));
    QTextEdit *mseErrorValue = new QTextEdit;
    mseErrorValue->setReadOnly(true);

    QVBoxLayout *valuesLayout = new QVBoxLayout;
    valuesLayout->addWidget(exactValueLabel);
    valuesLayout->addWidget(exactValue);
    valuesLayout->addWidget(approxValueLabel);
    valuesLayout->addWidget(approxValue);
    valuesGroupBox->setLayout(valuesLayout);
    valuesGroupBox->setFixedSize(320,140);
    valuesGroupBox->move(50,450);

    QVBoxLayout *errorLayout = new QVBoxLayout;
    errorLayout->addWidget(rErrorValueLabel);
    errorLayout->addWidget(rErrorValue);
    errorLayout->addWidget(mseValueLabel);
    errorLayout->addWidget(mseErrorValue);
    errorGroupBox->setLayout(errorLayout);
    errorGroupBox->setFixedSize(320,140);
    errorGroupBox->move(470,450);
}

void Window::createInputFieldsForParameters()
{
    QGroupBox *groupBox = new QGroupBox(tr("Parameters"), this);

    QLabel *maxIterLabel = new QLabel(tr("Max number of iteration (between %1 and %2)").arg(1).arg(10000));
    QSpinBox *maxIterSpinBox = new QSpinBox;
    maxIterSpinBox->setRange(1, 10000);
    maxIterSpinBox->setSingleStep(10);
    maxIterSpinBox->setValue(1000);

    QLabel *dLabel = new QLabel(tr("Space dimention d (between %1 and %2)").arg(1).arg(10));
    QSpinBox *dSpinBox = new QSpinBox;
    dSpinBox->setRange(1, 10);
    dSpinBox->setSingleStep(1);
    dSpinBox->setValue(3);

    QLabel *nLabel = new QLabel(tr("Number of function to interpolate"));
    QSpinBox *nSpinBox = new QSpinBox;
    nSpinBox->setRange(1, 100);
    nSpinBox->setSingleStep(1);
    nSpinBox->setValue(1);

    QVBoxLayout *parametersLayout = new QVBoxLayout;
    parametersLayout->addWidget(maxIterLabel);
    parametersLayout->addWidget(maxIterSpinBox);
    parametersLayout->addWidget(dLabel);
    parametersLayout->addWidget(dSpinBox);
    parametersLayout->addWidget(nLabel);
    parametersLayout->addWidget(nSpinBox);
    groupBox->setLayout(parametersLayout);

    groupBox->setFixedSize(320,200);
    groupBox->move(50,200);
}

void Window::createInputFieldForPointCoordinates()
{

    QGroupBox *groupBox = new QGroupBox(tr("Point coordinates X = (x,y,z)"), this);

    QLabel *xLabel = new QLabel(tr("x between %1 and %2:").arg(-1.0).arg(1.0));
    QDoubleSpinBox *xSpinBox = new QDoubleSpinBox;
    xSpinBox->setRange(-1.0, 1.0);
    xSpinBox->setSingleStep(0.1);
    xSpinBox->setValue(0);

    QLabel *yLabel = new QLabel(tr("y between %1 and %2:").arg(-1.0).arg(1.0));
    QDoubleSpinBox *ySpinBox = new QDoubleSpinBox;
    ySpinBox->setRange(-1.0, 1.0);
    ySpinBox->setSingleStep(0.1);
    ySpinBox->setValue(0);

    QLabel *zLabel = new QLabel(tr("z between %1 and %2:").arg(-1.0).arg(1.0));
    QDoubleSpinBox *zSpinBox = new QDoubleSpinBox;
    zSpinBox->setRange(-1.0, 1.0);
    zSpinBox->setSingleStep(0.1);
    zSpinBox->setValue(0);

    QVBoxLayout *coordinatesLayout = new QVBoxLayout;
    coordinatesLayout->addWidget(xLabel);
    coordinatesLayout->addWidget(xSpinBox);
    coordinatesLayout->addWidget(yLabel);
    coordinatesLayout->addWidget(ySpinBox);
    coordinatesLayout->addWidget(zLabel);
    coordinatesLayout->addWidget(zSpinBox);
    groupBox->setLayout(coordinatesLayout);

    groupBox->setFixedSize(320,200);
    groupBox->move(470,200);
}

void Window::createExclusiveGroupForFunctions()
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
}

void Window::createExclusiveGroupForMethods()
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

    groupBox->setFixedSize(420,100);
    groupBox->move(420,50);
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
