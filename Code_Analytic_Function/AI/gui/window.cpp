#include <QtWidgets>

#include "window.h"

std::string Window::projectPath = "/home/sialar/Stage/LaboJ_LLions/Code/";

Window::Window() : f(1), m0(0), m1(0), m2(0), m3(0), m4(0), d(3), n(1),
        maxIteration(1000), x(0), y(0), z(0), t(0), s(0), QWidget()
{
    int l = 1200, h = 800;
    int offset = 30;
    int lf = (l-4*offset)/3, hf = (h-3*offset)/2;
    int lr = lf, hr = hf;
    int lb = lf, hb = (h-8*offset)/2;
    int lp = lf, hp = hb;
    int lm = (lf-offset)/2, hm = h-10*offset;
    int lc = lm, hc = hm;
    setFixedSize(l, h);
    createExclusiveGroupForFunctions(lf,hf,offset,offset);
    createOutputFieldForResults(lr,hr,offset,2*offset+hf);
    createControlButtons(lb,hb,2*offset+lf,2*offset);
    createInputFieldsForParameters(lp,hp,2*offset+lr,h-2*offset-hp);
    createInputFieldForMethod(lm,hm,l-2*offset-2*lm,5*offset);
    createInputFieldForPointCoordinates(lc,hc,l-offset-lc,5*offset);
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

void Window::createControlButtons(int a, int b, int c, int d)
{
    QGroupBox *buttonsGroupBox = new QGroupBox(tr(""), this);

    QPushButton *randomButton = new QPushButton("Random point",this);
    QPushButton *startButton = new QPushButton("Interpolate",this);
    QPushButton *plotButton = new QPushButton("Plot results (1d) ",this);
    QPushButton *plotPath = new QPushButton("Plot interpolation points (2d ou 3d)",this);
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

    buttonsGroupBox->setFixedSize(a,b);
    buttonsGroupBox->move(c,d);
}

void Window::createOutputFieldForResults(int a, int b, int c, int d)
{
    QGroupBox *valuesGroupBox = new QGroupBox(tr("Results"), this);

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
    valuesLayout->addWidget(rErrorValueLabel);
    valuesLayout->addWidget(rErrorValue);
    valuesLayout->addWidget(mseValueLabel);
    valuesLayout->addWidget(mseErrorValue);
    valuesLayout->addWidget(timeLabel);
    valuesLayout->addWidget(execTimeValue);
    valuesGroupBox->setLayout(valuesLayout);

    valuesGroupBox->setFixedSize(a,b);
    valuesGroupBox->move(c,d);
}

void Window::createInputFieldsForParameters(int a, int b, int c, int d)
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

    groupBox->setFixedSize(a,b);
    groupBox->move(c,d);
}

void Window::createInputFieldForMethod(int a, int b, int c, int d)
{

    QGroupBox *groupBox = new QGroupBox(tr("Methods"), this);
    QVBoxLayout *coordinatesLayout = new QVBoxLayout;

    mxSpinBox = new QSpinBox;
    mxSpinBox->setRange(0, 2);
    mxSpinBox->setSingleStep(1);
    mxSpinBox->setValue(0);

    mySpinBox = new QSpinBox;
    mySpinBox->setRange(0, 2);
    mySpinBox->setSingleStep(1);
    mySpinBox->setValue(0);

    mzSpinBox = new QSpinBox;
    mzSpinBox->setRange(0, 2);
    mzSpinBox->setSingleStep(1);
    mzSpinBox->setValue(0);

    mtSpinBox = new QSpinBox;
    mtSpinBox->setRange(0, 2);
    mtSpinBox->setSingleStep(1);
    mtSpinBox->setValue(0);

    msSpinBox = new QSpinBox;
    msSpinBox->setRange(0, 2);
    msSpinBox->setSingleStep(1);
    msSpinBox->setValue(0);

    connect(mxSpinBox,SIGNAL(valueChanged(int)),this,SLOT(mx(int)));
    connect(mySpinBox,SIGNAL(valueChanged(int)),this,SLOT(my(int)));
    connect(mzSpinBox,SIGNAL(valueChanged(int)),this,SLOT(mz(int)));
    connect(mtSpinBox,SIGNAL(valueChanged(int)),this,SLOT(mt(int)));
    connect(msSpinBox,SIGNAL(valueChanged(int)),this,SLOT(ms(int)));

    coordinatesLayout->addWidget(mxSpinBox);
    coordinatesLayout->addWidget(mySpinBox);
    coordinatesLayout->addWidget(mzSpinBox);
    coordinatesLayout->addWidget(mtSpinBox);
    coordinatesLayout->addWidget(msSpinBox);
    groupBox->setLayout(coordinatesLayout);

    groupBox->setFixedSize(a,b);
    groupBox->move(c,d);
}

void Window::createInputFieldForPointCoordinates(int a, int b, int c, int d)
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

    groupBox->setFixedSize(a,b);
    groupBox->move(c,d);
}

void Window::createExclusiveGroupForFunctions(int a, int b, int c, int d)
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

    groupBox->setFixedSize(a,b);
    groupBox->move(c,d);
}

void Window::startInterpolation()
{
    string cmd = projectPath + "bin/TestX " + to_string(d) + " " + to_string(n) + " 100 " + \
            " " + to_string(maxIteration) + " " + to_string(f) + " " + to_string(x) + " " + \
            to_string(y) + " " + to_string(z) + " " + to_string(t) + " " + to_string(s) + \
            " " + to_string(m0) + " " + to_string(m1) + " " + to_string(m2) + \
            " " + to_string(m3) + " " + to_string(m4);
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
    string cmd = projectPath + "exec.sh PLOT " + to_string(m0) + " " + to_string(n) + " " + to_string(f);
    system(cmd.c_str());
    std::ifstream file(projectPath + "data/res.txt", ios::in);
    getline(file,relative_error);
    getline(file,mse_error);
    getline(file,exec_time);
    updateGUI();
}

void Window::plotInterpolationPoints()
{
    string cmd = projectPath + "exec.sh PATH " + to_string(d) + " 1 10 " + to_string(m0) + " " + \
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
