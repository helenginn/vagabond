//
//  Options.h
//  vagabond
//
//  Created by Helen Ginn on 22/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Options__
#define __vagabond__Options__

#include <stdio.h>
#include <vector>
#include <string>
#include "shared_ptrs.h"
#include "Crystal.h"
#include "Notifiable.h"

typedef enum
{
    ModelFilePDB,
    ModelFileVagabond,
} ModelFile;

class Options
{
public:
    Options(int argc, const char **argv);
    void run();

    static void setRuntimeOptions(OptionsPtr pointer)
    {
        Options::options = pointer;
        pointer->setManual(true);
    }

    static OptionsPtr getRuntimeOptions()
    {
        return options;
    }

    void setNotify(Notifiable *notifiable)
    {
        _notify = notifiable;
    }

    Notifiable *getNotify()
    {
        return _notify;
    }

    size_t crystalCount()
    {
        return crystals.size();
    }

    CrystalPtr getCrystal(int i)
    {
        return crystals[i];
    }

    CrystalPtr getActiveCrystal()
    {
        if (!crystals.size()) return CrystalPtr();
        return crystals[0];
    }

    static double getKick()
    {
        return _kick;
    }

    static double getDampen()
    {
        return _dampen;
    }

    static double minRes()
    {
        return _minRes;
    }

    static int enableTests()
    {
        return _enableTests;
    }

    static double getBStart()
    {
        return _bStart;
    }

    static double getBMult()
    {
        return _bMult;
    }
    
    static void setBMult(double bMult)
    {
        _bMult = bMult;
    }
    
    void setManual(bool manual)
    {
        _manual = manual;
    }

    void statusMessage(std::string message);
    void agreementSummary();
    void refineAll(RefinementType type, int numCycles, int *count = NULL,
            bool keepGoing = false);
    void superimposeAll(CrystalPtr crystal = CrystalPtr());
    void applyBMultiplier();
    void openModel(std::string pdbName);
    void openMTZ(std::string mtzName);
    void recalculateFFT();
private:
    static OptionsPtr options;
    Notifiable *_notify;
    void notifyGUI(bool enable);

    void parse();
    void displayHelp();
    void outputCrystalInfo();
    void refinementCycle(MoleculePtr molecule, int *count,
                             RefinementType type);
    bool parseJoke(std::string arg);

    std::vector<std::string> arguments;

    std::vector<ObjectPtr> objects;
    std::vector<CrystalPtr> crystals;
    std::vector<DatasetPtr> datasets;
    std::vector<DiffractionPtr> diffractions;

    int _numCycles;
    int _globalCount;
    bool _tie;
    bool _manual;

    static double _kick;
    static double _dampen;
    static double _bMult;
    static int _enableTests;
    static double _bStart;
    std::string _outputDir;
        static double _minRes;
};

#endif /* defined(__vagabond__Options__) */
