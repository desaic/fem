#include "Element.hpp"
#include "ElementCoarse.hpp"
#include "ElementRegGrid.hpp"

std::vector<std::vector<int> > makevifine();
std::vector<Vector3f> makeXfine();

std::vector<Vector3f> ElementCoarse::Xfine = makeXfine();
std::vector<std::vector<int> > ElementCoarse::vifine = makevifine();

extern int oneEleV [8][3];

std::vector<ElementCoarse *> coarsen(const ElementRegGrid & grid)
{
  int nx = grid.nx/2;
  int ny = grid.ny/2;
  int nz = grid.nz/2;
  int idx = 0;
  std::vector<ElementCoarse *> elts(nx*ny*nz,0);
  for(int ii = 0; ii<nx; ii++){
    for(int jj = 0; jj<ny; jj++){
      for(int kk = 0; kk<nz; kk++){
        ElementCoarse * ec = new ElementCoarse();
        for(int ll = 0; ll<8; ll++){
          ec->fineEle[ll] = grid.GetEleInd(2*ii+oneEleV[ll][0], 
            2*jj+oneEleV[ll][1], 2*kk+oneEleV[ll][2]);
        }
        elts[idx] = ec;
        idx++;
      }
    }
  }

  return elts;
}

Vector3f ElementCoarse::getDisp(Vector3f p)
{
  Vector3f u;

  return u;
}

ElementCoarse::ElementCoarse():fineEle(8,0)
{}

std::vector<Vector3f> makeXfine()
{
  //27 fine vertices in a coarse element.
  std::vector<Vector3f> Xfine(27);
  int idx = 0;
  for(int ii = -1; ii<2; ii++){
    for(int jj = -1; jj<2; jj++){
      for(int kk = -1; kk<2; kk++){
        Xfine[idx] = Vector3f((float)ii, (float)jj, (float)kk);
        idx++;
      }
    }
  }
  return Xfine;
}

std::vector<std::vector<int> > makevifine()
{
  //8 fine elements in a coarse element.
  std::vector<std::vector<int> > vifine(8);
  ElementRegGrid grid(2,2,2);
  for(unsigned int ii = 0; ii<grid.e.size();ii++){
    Element * e = grid.e[ii];
    vifine[ii].resize(e->nV());
    for(int jj = 0; jj<e->nV(); jj++){
      vifine[ii][jj] = e->at(jj);
    }
  }
  return vifine;
}