#include "memoryHandling.hh"

static MemoryList theList = 0;

PetscErrorCode registerPointerForDeallocation(void* thePointer)
{

  PetscErrorCode ierr;
  MemoryList temp;

  if (!theList){
    ierr = PetscMalloc1(1,&theList); CHKERRQ(ierr);
    theList->pointer = thePointer;
    theList->next = 0;
    ierr = PetscRegisterFinalize(deallocateAllRegisteredPointers); CHKERRQ(ierr);
  } else {
    // Traverse the list
    temp = theList;
    while (temp->next){
      temp = temp->next;
    }
    ierr = PetscMalloc1(1,&temp->next); CHKERRQ(ierr);
    temp = temp->next;
    temp->pointer = thePointer;
    temp->next = 0;
  }

  return 0;

}

PetscErrorCode deallocateAllRegisteredPointers()
{

  PetscErrorCode ierr;
  MemoryList temp;

  while (theList){
    ierr = PetscFree(theList->pointer); CHKERRQ(ierr);
    temp = theList;
    theList = theList->next;
    ierr = PetscFree(temp); CHKERRQ(ierr);
  }
  
  return 0;
  
}
