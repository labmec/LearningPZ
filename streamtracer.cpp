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
#include "TPZHDivApproxCreator.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "TPZMultiphysicsCompMesh.h"

enum EMatid {ENone,EDomain,EBC};
const int global_nthread = 16;

TPZGeoMesh* CreateGMesh(int ndiv, int dim);
TPZGeoMesh* ReadMeshFromGmsh(std::string file_name);
void CreateBCs(TPZGeoMesh* gmesh);
TPZCompMesh* CreateH1CMesh(TPZGeoMesh* gmesh, const int pord, TElasticity3DAnalytic *elas);

TPZCompMesh* CreateHdivCMesh(TPZGeoMesh* gmesh, const int pord);
void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh);
void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh);

void StreamTrace(TPZCompMesh* cmesh);

int main() {
    
    std::cout << "--------- Starting simulation ---------" << std::endl;
    const int pord = 1;
    int ndiv = 20;
    const int dim = 2;
    TPZGeoMesh* gmesh = CreateGMesh(ndiv,dim);
    std::ofstream out("gmesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);

    TPZCompMesh* cmeshHdiv = CreateHdivCMesh(gmesh,pord);

    TPZLinearAnalysis an(cmeshHdiv);
    SolveProblemDirect(an,cmeshHdiv);
    std::cout << "--------- PostProcess ---------" << std::endl;
    PrintResults(an,cmeshHdiv);
    
    StreamTrace(cmeshHdiv);

    delete cmeshHdiv;
    delete gmesh;
        
    std::cout << "--------- Simulation finished ---------" << std::endl;
}

TPZGeoMesh* CreateGMesh(int ndiv, int dim) {

    TPZGeoMesh* gmesh = new TPZGeoMesh;
    
    MMeshType meshType = MMeshType::ENoType;
    if(dim == 3)
        meshType = MMeshType::EHexahedral;
    if(dim == 2)
        meshType = MMeshType::EQuadrilateral;
    else
        DebugStop();        
   
    TPZManVector<REAL,3> minX = {-1,-1,-1};
    TPZManVector<REAL,3> maxX = {1,1,1};
    int nMats = 2*dim+1;
    
    constexpr bool createBoundEls{true};
    TPZVec<int> matIds(nMats,EBC);
    matIds[0] = EDomain;
    
    TPZManVector<int,3> ndivvec = {ndiv,ndiv,ndiv};
    if(dim == 2) ndivvec = {ndiv,ndiv};
    gmesh = TPZGeoMeshTools::CreateGeoMeshOnGrid(dim, minX, maxX, matIds, ndivvec, meshType, createBoundEls);

    return gmesh;
}


TPZCompMesh* CreateH1CMesh(TPZGeoMesh* gmesh, const int pord, TElasticity3DAnalytic *elas) {

    TPZCompMesh* cmesh = new TPZCompMesh(gmesh);
    const int dim = gmesh->Dimension();
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(pord);
    cmesh->SetAllCreateFunctionsContinuous();
    
    const STATE E = elas->fE, nu = elas->fPoisson;
    TPZManVector<STATE> force = {0,0,0};
    TPZElasticity3D *mat = new TPZElasticity3D(EDomain, E, nu, force, 0., 0., 0.);
    mat->SetExactSol(elas->ExactSolution(), 2);
    mat->SetForcingFunction(elas->ForceFunc(), 4);
    cmesh->InsertMaterialObject(mat);
    
    TPZFMatrix<STATE> val1(3,3,0.);
    TPZManVector<STATE> val2(3,0.);
    
    const int diri = 0, neu = 1, mixed = 2, normaltrac = 4;
    auto* BCCond0 = mat->CreateBC(mat, EBC, diri, val1, val2);
    BCCond0->SetForcingFunctionBC(elas->ExactSolution(), 4);
    cmesh->InsertMaterialObject(BCCond0);
    
    cmesh->AutoBuild();
    
    return cmesh;

}

void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh)

{
    TPZSkylineStructMatrix<STATE> matskl(cmesh);
    matskl.SetNumThreads(global_nthread);
    an.SetStructuralMatrix(matskl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    std::cout << "--------- Assemble ---------" << std::endl;
    TPZSimpleTimer time_ass;
    an.Assemble();
    std::cout << "Total time = " << time_ass.ReturnTimeDouble()/1000. << " s" << std::endl;
    std::cout << "--------- Solve ---------" << std::endl;
    TPZSimpleTimer time_sol;
    an.Solve();
    std::cout << "Total time = " << time_sol.ReturnTimeDouble()/1000. << " s" << std::endl;
    
    return;
}

void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh)

{
    std::cout << "--------- Post Process ---------" << std::endl;
    TPZSimpleTimer postProc("Post processing time");
    const std::string plotfile = "postprocess";
    constexpr int vtkRes{0};
    TPZVec<std::string> fields = {
        "Pressure",
        "Flux"
    };
    auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
    vtk.SetNThreads(global_nthread);
    vtk.Do();
    std::cout << "Total time = " << postProc.ReturnTimeDouble()/1000. << " s" << std::endl;
    
    return;
}

TPZCompMesh* CreateHdivCMesh(TPZGeoMesh* gmesh, const int pord) {

    TPZHDivApproxCreator hdivcreator(gmesh);
    hdivcreator.HdivFamily() = HDivFamily::EHDivStandard;
    hdivcreator.ProbType() = ProblemType::EDarcy;
    hdivcreator.IsRigidBodySpaces() = false;
    hdivcreator.SetDefaultOrder(pord);
    hdivcreator.SetExtraInternalOrder(0);
    hdivcreator.SetShouldCondense(false);
    hdivcreator.HybridType() = HybridizationType::ENone;

    const int dim = gmesh->Dimension();
    
    TPZMixedDarcyFlow *mat = new TPZMixedDarcyFlow(EDomain, dim);
    mat->SetConstantPermeability(1.); 

    auto forcingfunc = [](const TPZVec<REAL>& x, TPZVec<STATE>&f) {
        f[0] = x[0]*x[0] + x[1]*x[1] + sin(x[0]) + exp(x[1]);
    };
    auto forcingfuncbc = [&forcingfunc](const TPZVec<REAL>& x, TPZVec<STATE>&f, TPZMatrix<STATE>& mat) {
        forcingfunc(x,f);
    };
    mat->SetForcingFunction(forcingfunc,3);
    hdivcreator.InsertMaterialObject(mat);

    // Boundary conditions

    TPZFMatrix<STATE> val1(3,3,0.);
    TPZManVector<STATE> val2(3,1.);
    const int diri = 0, neu = 1;
    auto* BCCond0 = mat->CreateBC(mat, EBC, diri, val1, val2);
    BCCond0->SetForcingFunctionBC(forcingfuncbc, 4);
    hdivcreator.InsertMaterialObject(BCCond0);

    TPZMultiphysicsCompMesh *cmesh = hdivcreator.CreateApproximationSpace();
    return cmesh;

}

void StreamTrace(TPZCompMesh* cmesh) {
    // Aqui implementar funcao que segue o caminho de uma particula no dominio.
    // Deve comecar de um certo x,y e terminar em outro x,y.
    // Note que devera calcular como seguir uma particula para cada elemento, identificando onde ela entrou e saiu.
    // Depois, identificar qual o elemento vizinho onde agora entrou, e repetir o processo.
}