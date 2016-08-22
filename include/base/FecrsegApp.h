#ifndef FECRSEGAPP_H
#define FECRSEGAPP_H

#include "MooseApp.h"

class FecrsegApp;

template<>
InputParameters validParams<FecrsegApp>();

class FecrsegApp : public MooseApp
{
public:
  FecrsegApp(InputParameters parameters);
  virtual ~FecrsegApp();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
};

#endif /* FECRSEGAPP_H */
