#include "confdefs.h"
#include "petsc.h"
#if PETSC_VERSION_MINOR==0
#include "petscfix.h"
#include "petscis.h"
#endif

#include <string.h>

extern "C" {
  void iscopyindices_(IS *is, PetscInt *iarray,PetscErrorCode *ierr);
  void petscobjectreference_(PetscObject obj, int *__ierr );
  void pcfieldsplitgetsubksp_(PC *pc, KSP subksps[], int *__ierr);
}

// our own version of IsGetIndices, that just does a copy
// to avoid silly fortran interface issues
void iscopyindices_(IS *is, PetscInt *iarray,PetscErrorCode *ierr){
#if PETSC_VERSION_MAJOR==3
  const PetscInt *iptr;
#else
  PetscInt *iptr;
#endif
  PetscInt n;
  
  *ierr = ISGetIndices(*is, &iptr); if (*ierr) return;
  *ierr = ISGetSize(*is, &n); if (*ierr) return;
  
  memcpy(iarray, iptr, sizeof(PetscInt)*n);
  
  *ierr = ISRestoreIndices(*is, &iptr);
}

// PetscObjectReference()
// missing from the fortran interface (copy of petscobjectderefence, which
// does exists, with obvious changes)

#ifdef PETSC_USE_POINTER_CONVERSION
#if defined(__cplusplus)
extern "C" { 
#endif 
extern void *PetscToPointer(void*);
extern int PetscFromPointer(void *);
extern void PetscRmPointer(void*);
#if defined(__cplusplus)
} 
#endif 

#else

#define PetscToPointer(a) (*(long *)(a))
#define PetscFromPointer(a) (long)(a)
#define PetscRmPointer(a)
#endif

void petscobjectreference_(PetscObject obj, int *__ierr ){
*__ierr = PetscObjectReference(
   (PetscObject)PetscToPointer((obj) ));
}

// PCFieldSplitGetSubKSP()
// Again, missing from petsc fortrarn interface (sigh)
// Also, work around clumsy allocated subksps array that is returned
// by assuming we know how many ksp there are
void pcfieldsplitgetsubksp_(PC *pc, KSP subksps[], int *__ierr){
  KSP *temp_subksps;
  PetscInt nsubksps;
  *__ierr = PCFieldSplitGetSubKSP(*pc, &nsubksps, &temp_subksps);
  for (int i=0; i<nsubksps; i++) {
    subksps[i]=temp_subksps[i];
  }
  free(temp_subksps);
}
