#ifndef GLOBALDEFINITIONS
#define GLOBALDEFINITIONS

struct IncompressibleFlowParameters {

  IncompressibleFlowParameters() : velocityChar(), lengthChar(), ReynoldsNumber() {}
  PetscScalar velocityChar;
  PetscScalar lengthChar;
  PetscScalar ReynoldsNumber;

};

struct ThermalFlowParameters {

  ThermalFlowParameters()
    : incPar(), deltaT(), PrandtlNumber(), RayleighNumber() {}
  IncompressibleFlowParameters incPar;
  PetscScalar deltaT;
  PetscScalar PrandtlNumber;
  PetscScalar RayleighNumber;

};

// This should be moved somewhere more sensible
class BoxRange {

public:

  class Iterator {
    friend class BoxRange;
  public:
    PetscInt operator*() const { return index; }
    const Iterator& operator++() { ++index; return *this; }
    Iterator operator++(PetscInt){ Iterator copy(*this); ++index; return copy; }

    PetscBool operator==(const Iterator& rhs) const
    {
      return (PetscBool) (index == rhs.index);
    }
    PetscBool operator!=(const Iterator& rhs) const
    {
      return (PetscBool) (index <= rhs.index);
      /* This is a bit dirty, but it
         makes the range loops work, and 
         those minimize the chance I have to
         debug loops, and that, in turn, makes
         me a happy camper. I'll live with the
         smudge, is what I'm saying.
      */
    }

  protected:
    Iterator(PetscInt _ind) : index(_ind) {}

  private:
    PetscInt index;
  };

  Iterator begin() const { return Iterator(beginId); }
  Iterator end() const { return Iterator(endId); }
  PetscInt getBeginId() const { return beginId; }
  PetscInt getEndId() const { return endId; }

  BoxRange(PetscInt _beg, PetscInt _end) : beginId(_beg), endId(_end) {}
  BoxRange() : beginId(), endId() {}

private:

  PetscInt beginId, endId;

};

#endif
