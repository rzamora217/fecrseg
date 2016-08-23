#include "FecrsegApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ActionFactory.h"
#include "MooseSyntax.h"

// Modules
#include "PhaseFieldApp.h"
#include "HeatConductionApp.h"
#include "SolidMechanicsApp.h"


// Kernels
#include "FeCrIdealSink.h"            // lanlFeCr
#include "CahnHilliardD.h"            // lanlFeCr
#include "CahnHilliardBaseD.h"        // lanlFeCr
#include "DiffusionAdd.h"             // lanlFeCr
#include "DiffusionAdd2.h"            // lanlFeCr
#include "FeCrDRecombRate.h"          // lanlFeCr
#include "CHRadSinkFeCr.h"            // lanlFeCr
#include "CHRadSinkFeVa.h"            // lanlFeCr
#include "CHRadSinkFeSI.h"            // lanlFeCr
#include "CHRadSinkFeCrNoPF.h"        // lanlFeCr
#include "CHRadSinkFeVaNoPF.h"        // lanlFeCr
#include "CHRadSinkFeSINoPF.h"        // lanlFeCr
#include "UISource.h"                 // lanlFeCr
#include "VaSource.h"                 // lanlFeCr

// Materials
#include "FeCrVaSI.h"                 // lanlFeCr
#include "FeCrVaSIbulk.h"             // lanlFeCr

template<>
InputParameters validParams<FecrsegApp>()
{
  InputParameters params = validParams<MooseApp>();

  params.set<bool>("use_legacy_uo_initialization") = false;
  params.set<bool>("use_legacy_uo_aux_computation") = false;
  params.set<bool>("use_legacy_output_syntax") = false;

  return params;
}

FecrsegApp::FecrsegApp(InputParameters parameters) :
    MooseApp(parameters)
{
  Moose::registerObjects(_factory);
  PhaseFieldApp::registerObjects(_factory);
  HeatConductionApp::registerObjects(_factory);
  SolidMechanicsApp::registerObjects(_factory);
  FecrsegApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  PhaseFieldApp::associateSyntax(_syntax, _action_factory);
  HeatConductionApp::associateSyntax(_syntax, _action_factory);
  SolidMechanicsApp::associateSyntax(_syntax, _action_factory);
  FecrsegApp::associateSyntax(_syntax, _action_factory);
}

FecrsegApp::~FecrsegApp()
{
}

// External entry point for dynamic application loading
extern "C" void FecrsegApp__registerApps() { FecrsegApp::registerApps(); }
void
FecrsegApp::registerApps()
{
  registerApp(FecrsegApp);
}

// External entry point for dynamic object registration
extern "C" void FecrsegApp__registerObjects(Factory & factory) { FecrsegApp::registerObjects(factory); }
void
FecrsegApp::registerObjects(Factory & factory)
{
   // Kernels
  registerKernel(FeCrIdealSink);       // lanlFeCr
  registerKernel(CahnHilliardD);       // lanlFeCr
  registerKernel(DiffusionAdd);        // lanlFeCr
  registerKernel(DiffusionAdd2);       // lanlFeCr
  registerKernel(FeCrDRecombRate);     // lanlFeCr
  registerKernel(CHRadSinkFeCr);       // lanlFeCr
  registerKernel(CHRadSinkFeVa);       // lanlFeCr
  registerKernel(CHRadSinkFeSI);       // lanlFeCr
  registerKernel(CHRadSinkFeCrNoPF);   // lanlFeCr
  registerKernel(CHRadSinkFeVaNoPF);   // lanlFeCr
  registerKernel(CHRadSinkFeSINoPF);   // lanlFeCr
  registerKernel(UISource);            // lanlFeCr
  registerKernel(VaSource);            // lanlFeCr

  // Materials
  registerMaterial(FeCrVaSI);          // lanlFeCr
  registerMaterial(FeCrVaSIbulk);      // lanlFeCr
}

// External entry point for dynamic syntax association
extern "C" void FecrsegApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory) { FecrsegApp::associateSyntax(syntax, action_factory); }
void
FecrsegApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
}
