#include <QApplication>

#include "window.h"
#include <time.h>

int main(int argc, char *argv[])
{
    system("cd /home/sialar/Stage/LaboJ_LLions/Code/ \nmake -j8");
    srand(time(NULL));

    srand(time(NULL));
    QApplication app(argc, argv);
    Window w;
    w.show();
    return app.exec();
}
