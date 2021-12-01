QT       += core gui concurrent charts

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++11

INCLUDEPATH += D:/CppLibs/eigen-3.4.0

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

TARGET = UniSimulator.exe

RC_ICONS = icon.ico


SOURCES += \
    Simulator.cpp \
    WidgetCodes.cpp \
    main.cpp \
    MainWindow.cpp \
    tests.cpp

HEADERS += \
    MainWindow.h \
    Simulator.h \
    defines.h \
    tests.h

FORMS += \
    MainWindow.ui

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
