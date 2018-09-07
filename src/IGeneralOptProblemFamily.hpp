#pragma once

#include "IGeneralOptProblem.hpp"

/// Базовый класс для семейств задач оптимизации
class IGeneralOptProblemFamily
{
protected:
  // указатели на задачи семейства, задачи должны быть сформированы в конструкторе класса
  // число задач в семействе = размеру вектора pOptProblems
  vector<IGeneralOptProblem*> pOptProblems;

  // Конструктор скрыт, чтобы нельзя было создать объект данного класса
  // Требуется объявить потомка, перенести конструктор в public и реализовать
  IGeneralOptProblemFamily() {}
public:
  int GetFamilySize() const;
  ~IGeneralOptProblemFamily();
};
