#pragma once

#include "IGeneralOptProblemFamily.hpp"
#include "IOptProblem.hpp"

class IOptProblemFamily : public IGeneralOptProblemFamily
{
protected:
  // Конструктор скрыт, чтобы нельзя было создать объект данного класса
  // Требуется объявить потомка, перенести конструктор в public и реализовать
  IOptProblemFamily() {}
public:
  IOptProblem * operator[](int index);
};