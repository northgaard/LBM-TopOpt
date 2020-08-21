#include "binaryIO.hh"

PetscErrorCode writeVectorToBinary(const std::string& outputFolder,
                                   const std::string& fileName,
                                   Vec theVec, MPI_Comm communicator)
{
  PetscErrorCode ierr;
  PetscViewer view;
  PetscFunctionBeginUser;
  std::string output;
  output = outputFolder + fileName;
  ierr = PetscViewerBinaryOpen(communicator,output.c_str(),FILE_MODE_WRITE,&view);
  CHKERRQ(ierr);
  ierr = VecView(theVec,view); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&view); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode loadVectorFromBinary(const std::string& fileFolder,
                                    const std::string& fileName,
                                    Vec theVec, MPI_Comm communicator)
{
  PetscErrorCode ierr;
  PetscViewer view;
  PetscFunctionBeginUser;
  std::string input;
  input = fileFolder + fileName;
  ierr = PetscViewerBinaryOpen(communicator,input.c_str(),FILE_MODE_READ,&view);
  CHKERRQ(ierr);
  ierr = VecLoad(theVec,view); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&view); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
