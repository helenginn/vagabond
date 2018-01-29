/****************************************************************************
** Meta object code from reading C++ file 'VagWindow.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.10.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "VagWindow.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'VagWindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.10.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_VagWindow_t {
    QByteArrayData data[11];
    char stringdata0[143];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_VagWindow_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_VagWindow_t qt_meta_stringdata_VagWindow = {
    {
QT_MOC_LITERAL(0, 0, 9), // "VagWindow"
QT_MOC_LITERAL(1, 10, 15), // "pushSuperimpose"
QT_MOC_LITERAL(2, 26, 0), // ""
QT_MOC_LITERAL(3, 27, 19), // "pushRefinePositions"
QT_MOC_LITERAL(4, 47, 21), // "pushRefineFlexibility"
QT_MOC_LITERAL(5, 69, 15), // "pushBMultiplier"
QT_MOC_LITERAL(6, 85, 16), // "pushExploreMcule"
QT_MOC_LITERAL(7, 102, 14), // "recalculateFFT"
QT_MOC_LITERAL(8, 117, 7), // "openPDB"
QT_MOC_LITERAL(9, 125, 7), // "openMTZ"
QT_MOC_LITERAL(10, 133, 9) // "setOutput"

    },
    "VagWindow\0pushSuperimpose\0\0"
    "pushRefinePositions\0pushRefineFlexibility\0"
    "pushBMultiplier\0pushExploreMcule\0"
    "recalculateFFT\0openPDB\0openMTZ\0setOutput"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_VagWindow[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       9,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   59,    2, 0x08 /* Private */,
       3,    0,   60,    2, 0x08 /* Private */,
       4,    0,   61,    2, 0x08 /* Private */,
       5,    0,   62,    2, 0x08 /* Private */,
       6,    0,   63,    2, 0x08 /* Private */,
       7,    0,   64,    2, 0x08 /* Private */,
       8,    0,   65,    2, 0x08 /* Private */,
       9,    0,   66,    2, 0x08 /* Private */,
      10,    0,   67,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void VagWindow::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        VagWindow *_t = static_cast<VagWindow *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->pushSuperimpose(); break;
        case 1: _t->pushRefinePositions(); break;
        case 2: _t->pushRefineFlexibility(); break;
        case 3: _t->pushBMultiplier(); break;
        case 4: _t->pushExploreMcule(); break;
        case 5: _t->recalculateFFT(); break;
        case 6: _t->openPDB(); break;
        case 7: _t->openMTZ(); break;
        case 8: _t->setOutput(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject VagWindow::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_VagWindow.data,
      qt_meta_data_VagWindow,  qt_static_metacall, nullptr, nullptr}
};


const QMetaObject *VagWindow::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *VagWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_VagWindow.stringdata0))
        return static_cast<void*>(this);
    if (!strcmp(_clname, "Notifiable"))
        return static_cast< Notifiable*>(this);
    return QMainWindow::qt_metacast(_clname);
}

int VagWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 9)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 9;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 9)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 9;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
