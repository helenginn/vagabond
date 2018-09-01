#include "MoleculeExplorer.h"
#include "SequenceView.h"
#include "../../libsrc/Monomer.h"
#include "MonomerExplorer.h"

MoleculeExplorer::MoleculeExplorer(QWidget *parent, MoleculePtr molecule)
{
	_monomerExplorer = NULL; 
	_molecule = molecule;

	this->resize(600, 400);
	QString title = "Explore molecule ";
	title += QString::fromStdString(molecule->getChainID());
	this->setWindowTitle(title);

	_sequenceView = new SequenceView(this, molecule);
	_sequenceView->setExplorer(this);
	_sequenceView->show();

	_scrollArea = new QScrollArea(this);
	_scrollArea->setBackgroundRole(QPalette::Dark);
	_scrollArea->setWidget(_sequenceView);
	_scrollArea->setGeometry(0, 0, 600, 48);
	_scrollArea->show();
}

void MoleculeExplorer::setGLKeeper(GLKeeper *keeper)
{
	_keeper = keeper;
}       

void MoleculeExplorer::displayMonomer(MonomerPtr monomer)
{
	delete _monomerExplorer;
	_monomerExplorer = new MonomerExplorer(this, monomer);
	_monomerExplorer->setGeometry(0, 50, 600, 380);
	_monomerExplorer->setKeeper(_keeper);
	_monomerExplorer->show();
}

void MoleculeExplorer::updateCorrelation()
{
	if (_monomerExplorer)
	{
		_monomerExplorer->updateCorrelation(true);
	}
}

MoleculeExplorer::~MoleculeExplorer()
{
	/*
	delete _sequenceView;
	_sequenceView = NULL;
	*/

	delete _monomerExplorer;
	_monomerExplorer = NULL;

	delete _scrollArea;
	_scrollArea = NULL;
}
