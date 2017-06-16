#include <QApplication>

#include "window.h"
#include <time.h>

int main(int argc, char *argv[])
{
    srand(time(NULL));
    QApplication app(argc, argv);
    Window w;
    w.show();
    return app.exec();
}
