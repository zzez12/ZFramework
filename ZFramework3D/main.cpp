#include "zframework3d.h"
#include <QtGui/QApplication>
#include <QFile>
#include <QTextStream>
#include <QSettings>
#include "ZConsoleColor.h"

QString g_logFileName;
QString g_iniFileName;

void myMessageOutput(QtMsgType type, const char *msg)
{
	QString text;
	QTime time = QTime::currentTime();
	switch (type)
	{
	case QtDebugMsg:
		text = QString("%1--Debug: %2").arg(time.toString()).arg(msg);
		break;
	case QtWarningMsg:
		text = QString("%1--Warning: %2").arg(time.toString()).arg(msg);
		break;
	case QtCriticalMsg:
		text = QString("%1--Critical: %2").arg(time.toString()).arg(msg);
		break;
	case QtFatalMsg:
		text = QString("%1--Fatal: %2").arg(time.toString()).arg(msg);
		abort();
	}
	QFile file(g_logFileName);
	file.open(QIODevice::WriteOnly | QIODevice::Append);
	QTextStream ts(&file);
	std::cout << ZConsoleTools::red << text.toStdString() << std::endl;
	std::cout << ZConsoleTools::white;
	ts<<text<<endl;
}

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	g_logFileName = QCoreApplication::applicationDirPath() + QString("/zframework3d.log");
	g_iniFileName = QCoreApplication::applicationDirPath() + QString("/ZFramework3D.ini");
	// register log file
	qInstallMsgHandler(myMessageOutput);
 	QSettings *settings = new QSettings(g_iniFileName, QSettings::IniFormat);
	if (!settings)
		qDebug("Setting not loaded!");
	ZFramework3D w(settings);
	w.show();
	return a.exec();
}
