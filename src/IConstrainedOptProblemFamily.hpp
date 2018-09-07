#pragma once
#include "IConstrainedOptProblem.hpp"
#include "IGeneralOptProblemFamily.hpp"

class IConstrainedOptProblemFamily : public IGeneralOptProblemFamily
{
protected:
  // Конструктор скрыт, чтобы нельзя было создать объект данного класса
  // Требуется объявить потомка, перенести конструктор в public и реализовать
  IConstrainedOptProblemFamily() {}
public:
  IConstrainedOptProblem* operator[](int index);
};
