#include <QtWidgets>

#include "window.h"

std::string Window::projectPath = "/home/sialar/Stage/LaboJ_LLions/Code/";

Window::Window() : f(1), method(0), d(3), n(1),
        maxIteration(1000), x(0), y(0), z(0), t(0), s(0), QWidget()
{
    setFixedSize(1170, 650);
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
    execTimeValue->setText(QString::fromStdString(exec_time));
}

void Window::updateCoordsField()
{
    if (d<5) sSpinBox->setReadOnly(true);
    else sSpinBox->setReadOnly(false);
    if (d<4) tSpinBox->setReadOnly(true);
    else tSpinBox->setReadOnly(false);
    if (d<3) zSpinBox->setReadOnly(true);
    else zSpinBox->setReadOnly(false);
    if (d<2) ySpinBox->setReadOnly(true);
    else ySpinBox->setReadOnly(false);


    if (d<5) s=0;
    sSpinBox->setValue(s);
    if (d<4) t=0;
    else tSpinBox->setValue(t);
    if (d<3) z=0;
    else zSpinBox->setValue(z);
    if (d<2) y=0;
    else ySpinBox->setValue(y);
    xSpinBox->setValue(x);
}

void Window::createControlButtons()
{
    QGroupBox *buttonsGroupBox = new QGroupBox(tr(""), this);

    QPushButton *randomButton = new QPushButton("Random point",this);
    QPushButton *startButton = new QPushButton("Interpolate",this);
    QPushButton *plotButton = new QPushButton("Plot results",this);
    QPushButton *plotPath = new QPushButton("Plot interpolation points",this);
    QPushButton *errorButton = new QPushButton("Plot interpolation error",this);

    QObject::connect(randomButton, SIGNAL(clicked()), this, SLOT(randomPoint()));
    QObject::connect(startButton, SIGNAL(clicked()), this, SLOT(startInterpolation()));
    QObject::connect(plotPath, SIGNAL(clicked()), this, SLOT(plotInterpolationPoints()));
    QObject::connect(plotButton, SIGNAL(clicked()), this, SLOT(plotResults()));
    QObject::connect(errorButton, SIGNAL(clicked()), this, SLOT(plotErrors()));

    QVBoxLayout *buttonsLayout = new QVBoxLayout;
    buttonsLayout->addWidget(randomButton);
    buttonsLayout->addWidget(startButton);
    buttonsLayout->addWidget(plotButton);
    buttonsLayout->addWidget(plotPath);
    buttonsLayout->addWidget(errorButton);
    buttonsGroupBox->setLayout(buttonsLayout);

    buttonsGroupBox->setFixedSize(200,300);
    buttonsGroupBox->move(890, 160);
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

    QLabel *timeLabel = new QLabel(tr("Exectution time"));
    execTimeValue = new QTextEdit;
    execTimeValue->setReadOnly(true);

    QVBoxLayout *valuesLayout = new QVBoxLayout;
    valuesLayout->addWidget(exactValueLabel);
    valuesLayout->addWidget(exactValue);
    valuesLayout->addWidget(approxValueLabel);
    valuesLayout->addWidget(approxValue);
    valuesGroupBox->setLayout(valuesLayout);
    valuesGroupBox->setFixedSize(320,130);
    valuesGroupBox->move(470,470);

    QVBoxLayout *errorLayout = new QVBoxLayout;
    errorLayout->addWidget(rErrorValueLabel);
    errorLayout->addWidget(rErrorValue);
    errorLayout->addWidget(mseValueLabel);
    errorLayout->addWidget(mseErrorValue);
    errorLayout->addWidget(timeLabel);
    errorLayout->addWidget(execTimeValue);
    errorGroupBox->setLayout(errorLayout);
    errorGroupBox->setFixedSize(320,180);
    errorGroupBox->move(50,440);
}

void Window::createInputFieldsForParameters()
{
    QGroupBox *groupBox = new QGroupBox(tr("Parameters"), this);

    QLabel *maxIterLabel = new QLabel(tr("Max number of iteration (between %1 and %2)").arg(1).arg(10000));
    QSpinBox *maxIterSpinBox = new QSpinBox;
    maxIterSpinBox->setRange(1, 10000);
    maxIterSpinBox->setSingleStep(100);
    maxIterSpinBox->setValue(1000);

    QLabel *dLabel = new QLabel(tr("Space dimention d (between %1 and %2)").arg(1).arg(5));
    QSpinBox *dSpinBox = new QSpinBox;
    dSpinBox->setRange(1, 5);
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
    groupBox->move(50,190);
}

void Window::createInputFieldForPointCoordinates()
{

    QGroupBox *groupBox = new QGroupBox(tr("Point coordinates X"), this);
    QVBoxLayout *coordinatesLayout = new QVBoxLayout;

    QLabel *xLabel = new QLabel(tr("x between %1 and %2:").arg(-1.0).arg(1.0));
    xSpinBox = new QDoubleSpinBox;
    xSpinBox->setRange(-1.0, 1.0);
    xSpinBox->setSingleStep(0.1);
    xSpinBox->setValue(0);

    QLabel *yLabel = new QLabel(tr("y between %1 and %2:").arg(-1.0).arg(1.0));
    ySpinBox = new QDoubleSpinBox;
    ySpinBox->setRange(-1.0, 1.0);
    ySpinBox->setSingleStep(0.1);
    ySpinBox->setValue(0);

    QLabel *zLabel = new QLabel(tr("z between %1 and %2:").arg(-1.0).arg(1.0));
    zSpinBox = new QDoubleSpinBox;
    zSpinBox->setRange(-1.0, 1.0);
    zSpinBox->setSingleStep(0.1);
    zSpinBox->setValue(0);

    QLabel *tLabel = new QLabel(tr("t between %1 and %2:").arg(-1.0).arg(1.0));
    tSpinBox = new QDoubleSpinBox;
    tSpinBox->setRange(-1.0, 1.0);
    tSpinBox->setSingleStep(0.1);
    tSpinBox->setValue(0);

    QLabel *sLabel = new QLabel(tr("s between %1 and %2:").arg(-1.0).arg(1.0));
    sSpinBox = new QDoubleSpinBox;
    sSpinBox->setRange(-1.0, 1.0);
    sSpinBox->setSingleStep(0.1);
    sSpinBox->setValue(0);

    updateCoordsField();

    connect(xSpinBox,SIGNAL(valueChanged(double)),this,SLOT(updateX(double)));
    connect(ySpinBox,SIGNAL(valueChanged(double)),this,SLOT(updateY(double)));
    connect(zSpinBox,SIGNAL(valueChanged(double)),this,SLOT(updateZ(double)));
    connect(tSpinBox,SIGNAL(valueChanged(double)),this,SLOT(updateT(double)));
    connect(tSpinBox,SIGNAL(valueChanged(double)),this,SLOT(updateS(double)));

    coordinatesLayout->addWidget(xLabel);
    coordinatesLayout->addWidget(xSpinBox);
    coordinatesLayout->addWidget(yLabel);
    coordinatesLayout->addWidget(ySpinBox);
    coordinatesLayout->addWidget(zLabel);
    coordinatesLayout->addWidget(zSpinBox);
    coordinatesLayout->addWidget(tLabel);
    coordinatesLayout->addWidget(tSpinBox);
    coordinatesLayout->addWidget(sLabel);
    coordinatesLayout->addWidget(sSpinBox);
    groupBox->setLayout(coordinatesLayout);

    groupBox->setFixedSize(320,300);
    groupBox->move(470,160);
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

void Window::plotResults()
{
    int n = (maxIteration > 200) ? 100 : maxIteration;
    string cmd = projectPath + "exec.sh PLOT " + to_string(method) + " " + to_string(n) + " " + to_string(f);
    system(cmd.c_str());
    std::ifstream file(projectPath + "data/res.txt", ios::in);
    getline(file,relative_error);
    getline(file,mse_error);
    getline(file,exec_time);
    updateGUI();
}

void Window::plotInterpolationPoints()
{
    string cmd = projectPath + "exec.sh PATH " + to_string(d) + " 1 10 " + to_string(method) + " " + \
            to_string(maxIteration) + " " + to_string(f);
    std::cout << cmd << std::endl;
    system(cmd.c_str());
    std::ifstream file(projectPath + "data/res.txt", ios::in);
    getline(file,relative_error);
    getline(file,mse_error);
    getline(file,exec_time);
    updateGUI();
}

void Window::plotErrors()
{
    string cmd = projectPath + "exec.sh ERROR " + to_string(d) + " " + to_string(n) + " 10 ALL " + \
            to_string(maxIteration) + " 0 " + to_string(f);
    std::cout << cmd << std::endl;
    system(cmd.c_str());
    std::ifstream file(projectPath + "data/res.txt", ios::in);
    getline(file,relative_error);
    getline(file,mse_error);
    getline(file,exec_time);
    updateGUI();
}

void Window::randomPoint()
{
    x = (rand()/(double)RAND_MAX)*2 - 1;
    y = (rand()/(double)RAND_MAX)*2 - 1;
    z = (rand()/(double)RAND_MAX)*2 - 1;
    t = (rand()/(double)RAND_MAX)*2 - 1;
    s = (rand()/(double)RAND_MAX)*2 - 1;
    updateCoordsField();
}

Window::~Window()
{

}
