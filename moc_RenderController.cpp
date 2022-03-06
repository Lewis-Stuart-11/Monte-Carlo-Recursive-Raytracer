/****************************************************************************
** Meta object code from reading C++ file 'RenderController.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.13.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include <memory>
#include "RenderController.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'RenderController.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.13.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_RenderController_t {
    QByteArrayData data[33];
    char stringdata0[585];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_RenderController_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_RenderController_t qt_meta_stringdata_RenderController = {
    {
QT_MOC_LITERAL(0, 0, 16), // "RenderController"
QT_MOC_LITERAL(1, 17, 21), // "objectRotationChanged"
QT_MOC_LITERAL(2, 39, 0), // ""
QT_MOC_LITERAL(3, 40, 20), // "lightRotationChanged"
QT_MOC_LITERAL(4, 61, 11), // "zoomChanged"
QT_MOC_LITERAL(5, 73, 5), // "value"
QT_MOC_LITERAL(6, 79, 17), // "xTranslateChanged"
QT_MOC_LITERAL(7, 97, 17), // "yTranslateChanged"
QT_MOC_LITERAL(8, 115, 23), // "mapUVWToRGBCheckChanged"
QT_MOC_LITERAL(9, 139, 5), // "state"
QT_MOC_LITERAL(10, 145, 23), // "useLightingCheckChanged"
QT_MOC_LITERAL(11, 169, 29), // "texturedRenderingCheckChanged"
QT_MOC_LITERAL(12, 199, 29), // "textureModulationCheckChanged"
QT_MOC_LITERAL(13, 229, 21), // "depthTestCheckChanged"
QT_MOC_LITERAL(14, 251, 20), // "showAxesCheckChanged"
QT_MOC_LITERAL(15, 272, 22), // "showObjectCheckChanged"
QT_MOC_LITERAL(16, 295, 24), // "centreObjectCheckChanged"
QT_MOC_LITERAL(17, 320, 23), // "scaleObjectCheckChanged"
QT_MOC_LITERAL(18, 344, 14), // "shadowsChanged"
QT_MOC_LITERAL(19, 359, 18), // "reflectionsChanged"
QT_MOC_LITERAL(20, 378, 17), // "scatteringChanged"
QT_MOC_LITERAL(21, 396, 17), // "numSamplesChanged"
QT_MOC_LITERAL(22, 414, 20), // "emissiveLightChanged"
QT_MOC_LITERAL(23, 435, 19), // "ambientLightChanged"
QT_MOC_LITERAL(24, 455, 19), // "diffuseLightChanged"
QT_MOC_LITERAL(25, 475, 20), // "specularLightChanged"
QT_MOC_LITERAL(26, 496, 23), // "specularExponentChanged"
QT_MOC_LITERAL(27, 520, 15), // "BeginScaledDrag"
QT_MOC_LITERAL(28, 536, 11), // "whichButton"
QT_MOC_LITERAL(29, 548, 1), // "x"
QT_MOC_LITERAL(30, 550, 1), // "y"
QT_MOC_LITERAL(31, 552, 18), // "ContinueScaledDrag"
QT_MOC_LITERAL(32, 571, 13) // "EndScaledDrag"

    },
    "RenderController\0objectRotationChanged\0"
    "\0lightRotationChanged\0zoomChanged\0"
    "value\0xTranslateChanged\0yTranslateChanged\0"
    "mapUVWToRGBCheckChanged\0state\0"
    "useLightingCheckChanged\0"
    "texturedRenderingCheckChanged\0"
    "textureModulationCheckChanged\0"
    "depthTestCheckChanged\0showAxesCheckChanged\0"
    "showObjectCheckChanged\0centreObjectCheckChanged\0"
    "scaleObjectCheckChanged\0shadowsChanged\0"
    "reflectionsChanged\0scatteringChanged\0"
    "numSamplesChanged\0emissiveLightChanged\0"
    "ambientLightChanged\0diffuseLightChanged\0"
    "specularLightChanged\0specularExponentChanged\0"
    "BeginScaledDrag\0whichButton\0x\0y\0"
    "ContinueScaledDrag\0EndScaledDrag"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_RenderController[] = {

 // content:
       8,       // revision
       0,       // classname
       0,    0, // classinfo
      26,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,  144,    2, 0x0a /* Public */,
       3,    0,  145,    2, 0x0a /* Public */,
       4,    1,  146,    2, 0x0a /* Public */,
       6,    1,  149,    2, 0x0a /* Public */,
       7,    1,  152,    2, 0x0a /* Public */,
       8,    1,  155,    2, 0x0a /* Public */,
      10,    1,  158,    2, 0x0a /* Public */,
      11,    1,  161,    2, 0x0a /* Public */,
      12,    1,  164,    2, 0x0a /* Public */,
      13,    1,  167,    2, 0x0a /* Public */,
      14,    1,  170,    2, 0x0a /* Public */,
      15,    1,  173,    2, 0x0a /* Public */,
      16,    1,  176,    2, 0x0a /* Public */,
      17,    1,  179,    2, 0x0a /* Public */,
      18,    1,  182,    2, 0x0a /* Public */,
      19,    1,  185,    2, 0x0a /* Public */,
      20,    1,  188,    2, 0x0a /* Public */,
      21,    1,  191,    2, 0x0a /* Public */,
      22,    1,  194,    2, 0x0a /* Public */,
      23,    1,  197,    2, 0x0a /* Public */,
      24,    1,  200,    2, 0x0a /* Public */,
      25,    1,  203,    2, 0x0a /* Public */,
      26,    1,  206,    2, 0x0a /* Public */,
      27,    3,  209,    2, 0x0a /* Public */,
      31,    2,  216,    2, 0x0a /* Public */,
      32,    2,  221,    2, 0x0a /* Public */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int,    5,
    QMetaType::Void, QMetaType::Int,    5,
    QMetaType::Void, QMetaType::Int,    5,
    QMetaType::Void, QMetaType::Int,    9,
    QMetaType::Void, QMetaType::Int,    9,
    QMetaType::Void, QMetaType::Int,    9,
    QMetaType::Void, QMetaType::Int,    9,
    QMetaType::Void, QMetaType::Int,    9,
    QMetaType::Void, QMetaType::Int,    9,
    QMetaType::Void, QMetaType::Int,    9,
    QMetaType::Void, QMetaType::Int,    9,
    QMetaType::Void, QMetaType::Int,    9,
    QMetaType::Void, QMetaType::Int,    9,
    QMetaType::Void, QMetaType::Int,    9,
    QMetaType::Void, QMetaType::Int,    9,
    QMetaType::Void, QMetaType::Int,    5,
    QMetaType::Void, QMetaType::Int,    5,
    QMetaType::Void, QMetaType::Int,    5,
    QMetaType::Void, QMetaType::Int,    5,
    QMetaType::Void, QMetaType::Int,    5,
    QMetaType::Void, QMetaType::Int,    5,
    QMetaType::Void, QMetaType::Int, QMetaType::Float, QMetaType::Float,   28,   29,   30,
    QMetaType::Void, QMetaType::Float, QMetaType::Float,   29,   30,
    QMetaType::Void, QMetaType::Float, QMetaType::Float,   29,   30,

       0        // eod
};

void RenderController::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        auto *_t = static_cast<RenderController *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->objectRotationChanged(); break;
        case 1: _t->lightRotationChanged(); break;
        case 2: _t->zoomChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 3: _t->xTranslateChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 4: _t->yTranslateChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 5: _t->mapUVWToRGBCheckChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 6: _t->useLightingCheckChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 7: _t->texturedRenderingCheckChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 8: _t->textureModulationCheckChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 9: _t->depthTestCheckChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 10: _t->showAxesCheckChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 11: _t->showObjectCheckChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 12: _t->centreObjectCheckChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 13: _t->scaleObjectCheckChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 14: _t->shadowsChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 15: _t->reflectionsChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 16: _t->scatteringChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 17: _t->numSamplesChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 18: _t->emissiveLightChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 19: _t->ambientLightChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 20: _t->diffuseLightChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 21: _t->specularLightChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 22: _t->specularExponentChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 23: _t->BeginScaledDrag((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< float(*)>(_a[2])),(*reinterpret_cast< float(*)>(_a[3]))); break;
        case 24: _t->ContinueScaledDrag((*reinterpret_cast< float(*)>(_a[1])),(*reinterpret_cast< float(*)>(_a[2]))); break;
        case 25: _t->EndScaledDrag((*reinterpret_cast< float(*)>(_a[1])),(*reinterpret_cast< float(*)>(_a[2]))); break;
        default: ;
        }
    }
}

QT_INIT_METAOBJECT const QMetaObject RenderController::staticMetaObject = { {
    &QObject::staticMetaObject,
    qt_meta_stringdata_RenderController.data,
    qt_meta_data_RenderController,
    qt_static_metacall,
    nullptr,
    nullptr
} };


const QMetaObject *RenderController::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *RenderController::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_RenderController.stringdata0))
        return static_cast<void*>(this);
    return QObject::qt_metacast(_clname);
}

int RenderController::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 26)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 26;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 26)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 26;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
