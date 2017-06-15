#include <QtWidgets>

#include "window.h"

std::string Window::projectPath = "/home/sialar/Stage/LaboJ_LLions/Code/";

Window::Window() : f(1), method(0), d(3), n(1),
        maxIteration(1000), x(0), y(0), z(0), QWidget()
{
    setFixedSize(1150, 650);
    createExclusiveGroupForFunctions();
    createExclusiveGroupForMethods();
    createInputFieldForPointCoordinates();
    createInputFieldsForParameters();
    createOutputFieldForResults();
    createControlButtons();
}

void Window::updateGUI()
{
    exactValue->setText(QString::fromStdString(real_f));
    approxValue->setText(QString::fromStdString(f_tilde));
    rErrorValue->setText(QString::fromStdString(relative_error));
    mseErrorValue->setText(QString::fromStdString(mse_error));
}

void Window::createControlButtons()
{
    QGroupBox *buttonsGroupBox = new QGroupBox(tr(""), this);

    QPushButton *startButton = new QPushButton("Interpolate",this);
    QPushButton *plotButton = new QPushButton("Plot results",this);
    QPushButton *plotPath = new QPushButton("Plot interpolation points",this);
    QPushButton *displayResults = new QPushButton("Display results",this);

    QObject::connect(startButton, SIGNAL(clicked()), this, SLOT(startInterpolation()));

    QVBoxLayout *buttonsLayout = new QVBoxLayout;
    buttonsLayout->addWidget(startButton);
    buttonsLayout->addWidget(plotButton);
    buttonsLayout->addWidget(plotPath);
    buttonsLayout->addWidget(displayResults);
    buttonsGroupBox->setLayout(buttonsLayout);

    buttonsGroupBox->setFixedSize(200,200);
    buttonsGroupBox->move(890, 200);
}

void Window::createOutputFieldForResults()
{
    QGroupBox *valuesGroupBox = new QGroupBox(tr("Results"), this);
    QGroupBox *errorGroupBox = new QGroupBox(tr("Interpolation Error"), this);

    QLabel *exactValueLabel = new QLabel(tr("Exact value of f(X)"));
    exactValue = new QTextEdit;
    exactValue->setReadOnly(true);

    QLabel *approxValueLabel = new QLabel(tr("Approx value of f(X)"));
    approxValue = new QTextEdit;
    approxValue->setReadOnly(true);

    QLabel *rErrorValueLabel = new QLabel(tr("Relative interpolation error"));
    rErrorValue = new QTextEdit;
    rErrorValue->setReadOnly(true);

    QLabel *mseValueLabel = new QLabel(tr("MSE interpolation error"));
    mseErrorValue = new QTextEdit;
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

    connect(maxIterSpinBox,SIGNAL(valueChanged(int)),this,SLOT(updateMaxIter(int)));
    connect(dSpinBox,SIGNAL(valueChanged(int)),this,SLOT(updateD(int)));
    connect(nSpinBox,SIGNAL(valueChanged(int)),this,SLOT(updateN(int)));


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

    connect(xSpinBox,SIGNAL(valueChanged(double)),this,SLOT(updateX(double)));
    connect(ySpinBox,SIGNAL(valueChanged(double)),this,SLOT(updateY(double)));
    connect(zSpinBox,SIGNAL(valueChanged(double)),this,SLOT(updateZ(double)));

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

    connect(radio1, SIGNAL(clicked()), this, SLOT(f1()));
    connect(radio2, SIGNAL(clicked()), this, SLOT(f2()));
    connect(radio3, SIGNAL(clicked()), this, SLOT(f3()));

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

    connect(radio1, SIGNAL(clicked()), this, SLOT(m0()));
    connect(radio2, SIGNAL(clicked()), this, SLOT(m1()));
    connect(radio3, SIGNAL(clicked()), this, SLOT(m2()));

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

void Window::startInterpolation()
{
    string cmd = projectPath + "bin/TestX " + to_string(d) + " " + to_string(n) + " 100 " + \
            to_string(method) + " " + to_string(maxIteration) + " " + to_string(x) + " " + \
            to_string(y) + " " + to_string(z) + " " + to_string(f);
    std::cout << cmd << std::endl;
    system(cmd.c_str());
    std::ifstream file(projectPath + "data/res.txt", ios::in);

    getline(file,real_f);
    getline(file,f_tilde);
    getline(file,relative_error);
    getline(file,mse_error);
    getline(file,exec_time);

    updateGUI();
}

Window::~Window()
{

}
