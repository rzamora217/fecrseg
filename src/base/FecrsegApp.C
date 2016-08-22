#include "FecrsegApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

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
  ModulesApp::registerObjects(_factory);
  FecrsegApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
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
}

// External entry point for dynamic syntax association
extern "C" void FecrsegApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory) { FecrsegApp::associateSyntax(syntax, action_factory); }
void
FecrsegApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
}
