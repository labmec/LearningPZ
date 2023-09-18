#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <iostream>
#include "pzlog.h"
#include "pzgmesh.h"
#include "TPZGenGrid2D.h"
#include "TPZVTKGeoMesh.h"
#include "pzcmesh.h"
#include <DarcyFlow/TPZMixedDarcyFlow.h>
#include <TPZNullMaterial.h>
#include <pzbuildmultiphysicsmesh.h>
#include <pzskylstrmatrix.h>
#include <pzskylstrmatrix.h>
#include <pzstepsolver.h>
#include <TPZLinearAnalysis.h>
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#include <TPZSimpleTimer.h>
#include "pzvisualmatrix.h"
#include "TPZSYSMPMatrix.h"
#include "TPZVTKGenerator.h"
#include <Elasticity/TPZElasticity3D.h>
#include "TPZAnalyticSolution.h"
#include "TPZGeoMeshTools.h"
#include <TPZGmshReader.h>
#include "tpzchangeel.h"
#include "TPZRefPatternDataBase.h"
#include "TPZRefPatternTools.h"

enum EMatid {ENone,EDomain,EBC};
// enum: é um tipo de variável que seleciona um intervalo de valores;
// EMatid:
const int global_nthread = 16;


TPZGeoMesh* CreateGMesh(int ndiv);
TPZGeoMesh* ReadMeshFromGmsh(std::string file_name);
void CreateBCs(TPZGeoMesh* gmesh);
TPZCompMesh* CreateH1CMesh(TPZGeoMesh* gmesh, const int pord, TElasticity3DAnalytic *elas);

void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh);
void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh);


int main() {
    
    std::cout << "--------- Starting simulation ---------" << std::endl;
    //Imprime uma linha falando que a simulação se iniciou
    const int pord = 1;
    //declara a constante pord, utilizada para
    int ndiv = 2;
    //declara a constante ndiv, utilizada para a função CreateGMesh que armazena a memória de TPZGeoMesh
    TPZGeoMesh* gmesh = CreateGMesh(ndiv);
        
    std::ofstream out("gmesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);

    TElasticity3DAnalytic *elas = new TElasticity3DAnalytic;
    elas->fE = 250.;//206.8150271873455;
    elas->fPoisson = 0.;
    elas->fProblemType = TElasticity3DAnalytic::EStretchx;
    TPZCompMesh* cmeshH1 = CreateH1CMesh(gmesh,pord,elas);
    
    TPZLinearAnalysis an(cmeshH1);
    SolveProblemDirect(an,cmeshH1);
    //declara a constante pord, utilizada paras ---------" << std::endl;
    PrintResults(an,cmeshH1);
    
    // deleting stuff
    delete cmeshH1;
    delete gmesh;
        
    std::cout << "--------- Simulation finished ---------" << std::endl;
}pord

TPZGeoMesh* CreateGMesh(int ndiv) {
    //A função CreateGMesh serve para criar um endereço da memória onde o TPZGeoMesh está armazenado, a partir de um inteiro ndiv.
    TPZGeoMesh* gmesh = new TPZGeoMesh;
    // o gmesh armazena a memória onde o TPZGeoMesh está armazenado, servindo como um ponteiro para este objeto.
    MMeshType meshType = MMeshType::EHexahedral;
    int dim = 3;
    //declara a constante dim, utilizada para
    TPZManVector<REAL,3> minX = {-1,-1,-1};
    TPZManVector<REAL,3> maxX = {1,1,1};
    int nMats = 2*dim+1;
    //declara a constante nMats, definida por (2*dim+1), utilizada para
    
    constexpr bool createBoundEls{true};
    TPZVec<int> matIds(nMats,EBC);
    matIds[0] = EDomain;
    
    TPZManVector<int,3> ndivvec = {ndiv,ndiv,ndiv};
    gmesh = TPZGeoMeshTools::CreateGeoMeshOnGrid(dim, minX, maxX,matIds, ndivvec, meshType,createBoundEls);
    
    return gmesh;
}


TPZCompMesh* CreateH1CMesh(TPZGeoMesh* gmesh, const int pord, TElasticity3DAnalytic *elas) {
    TPZCompMesh* cmesh = new TPZCompMesh(gmesh);
    const int dim = gmesh->Dimension();
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(pord);
    cmesh->SetAllCreateFunctionsContinuous();
    
    // Domain elas mat
    const STATE E = elas->fE, nu = elas->fPoisson;
    TPZManVector<STATE> force = {0,0,0};
    TPZElasticity3D *mat = new TPZElasticity3D(EDomain, E, nu, force, 0., 0., 0.);
    mat->SetExactSol(elas->ExactSolution(), 2);
    mat->SetForcingFunction(elas->ForceFunc(), 4);
    cmesh->InsertMaterialObject(mat);
    
    // BC null mat
    TPZFMatrix<STATE> val1(3,3,0.);
    TPZManVector<STATE> val2(3,0.);
    
    const int diri = 0, neu = 1, mixed = 2, normaltrac = 4;
    auto* BCCond0 = mat->CreateBC(mat, EBC, diri, val1, val2);
    BCCond0->SetForcingFunctionBC(elas->ExactSolution(), 4);
    cmesh->InsertMaterialObject(BCCond0);
    
    // Constructs mesh
    cmesh->AutoBuild();
    
    return cmesh;
}

void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh)
{

    TPZSkylineStructMatrix<STATE> matskl(cmesh);
    // TPZSSpStructMatrix<STATE> matskl(cmesh);
    matskl.SetNumThreads(global_nthread);
    an.SetStructuralMatrix(matskl);
    
    ///Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);//ELU //ECholesky // ELDLt
    an.SetSolver(step);
    
    //assembles the system
    std::cout << "--------- Assemble ---------" << std::endl;
    TPZSimpleTimer time_ass;
    an.Assemble();
    std::cout << "Total time = " << time_ass.ReturnTimeDouble()/1000. << " s" << std::endl;
    
    ///solves the system
    std::cout << "--------- Solve ---------" << std::endl;
    TPZSimpleTimer time_sol;
    an.Solve();
    std::cout << "Total time = " << time_sol.ReturnTimeDouble()/1000. << " s" << std::endl;
    
    return;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh)
{
 
    std::cout << "--------- Post Process ---------" << std::endl;
    TPZSimpleTimer postProc("Post processing time");
    const std::string plotfile = "postprocess";
    constexpr int vtkRes{0};
    
    TPZVec<std::string> fields = {
        // "ExactDisplacement",
        // "ExactStress",
        "Displacement",
        "Stress",
        "Strain",
    };
    auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
    vtk.SetNThreads(global_nthread);
    vtk.Do();
    std::cout << "Total time = " << postProc.ReturnTimeDouble()/1000. << " s" << std::endl;
    
    return;
}
