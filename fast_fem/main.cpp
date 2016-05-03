//#####################################################################
// Elastic Hexahedron FEM
// Bo Zhu, MIT, 07/15
//#####################################################################
#include <iostream>
#include "PARSE_ARGS.h"
#include "Common.h"
#include "File.h"
#include "TopoOptimizationHelper.h"
#include "ElasticHexFEM.h"
#include "ElasticTetFEM.h"
#include "ElasticHexTetFEM.h"
#include "Homogenization.h"
#include "MicrostructureOpt.h"

#ifndef __Main_cpp__
#define __Main_cpp__

static const int d=3;
Type_Define_VectorD_Outside_Class(d);Type_Define_VectorDi_Outside_Class(d);Type_Define_MatrixD_Outside_Class(d); 

void TestHexFEM(PARSE_ARGS& parse_args)
{
	std::string output_dir=parse_args.Get_String_Value("-o");output_dir+="/0";
    const int s=parse_args.Get_Integer_Value("-s");
    const int solver_type=parse_args.Get_Integer_Value("-solver");
	const int processor_type=parse_args.Get_Integer_Value("-proc");
    const bool write_residual=parse_args.Get_Option_Value("-wr");

	T length=(T)1;VectorDi cell_counts=VectorDi::Ones()*s;
	if(d>1)cell_counts[1]=cell_counts[0]/2;if(d>2)cell_counts[2]=cell_counts[0]/2;
    Lattice<d> lattice(cell_counts,(T)1/(T)cell_counts[0]);
	ElasticHexFEM<d> elas(lattice);
	elas.solver_type=ElasticHexFEM<d>::SolverType(solver_type);
	elas.processor_type=ElasticHexFEM<d>::ProcessorType(processor_type);
	elas.des_iter_num=parse_args.Get_Integer_Value("-n_d");
	elas.asc_iter_num=parse_args.Get_Integer_Value("-n_a");
	elas.v_cycle_iter_num=parse_args.Get_Integer_Value("-n_c");
	elas.krylov_iter_num=parse_args.Get_Integer_Value("-krylov_iter_num");
	int material_pattern=parse_args.Get_Integer_Value("-mt");
	elas.use_heat=false;
	elas.use_implicit_psi_P=false;

	////reinitialize materials
	elas.materials.clear();elas.materials.push_back(ElasticMaterial((T)1,(T).3));
	for(LatticeIterator<d> iter(lattice,LatticeIterator<d>::CELL);iter.Valid();iter.Next()){elas.material_id(iter.Coord())=0;}

	//////reinitialize materials
	//for(LatticeIterator<d> iter(lattice,LatticeIterator<d>::CELL);iter.Valid();iter.Next()){
	//	VectorDi cell=iter.Coord();
	//	int c=0;for(int i=0;i<d;i++){
	//		if(cell[i]==0||cell[i]==lattice.cell_counts[i]-1)c++;
	//	}
	//	elas.material_id(cell)=(c>=d-1?0:-1);

	//	//switch(material_pattern){
	//	//case 0:{elas.material_id(cell)=0;}break;
	//	//case 1:{elas.material_id(cell)=(cell[0]<cell_counts[0]/2?0:1);}break;
	//	//case 2:{int s=16;elas.material_id(cell)=((cell[0]/s)%2==0?0:1);}break;
	//	//case 3:{int s=5;elas.material_id(cell)=((cell[1]/s)%2==0?0:1);}break;
	//	//case 4:{
	//	//	VectorDi offset=lattice.cell_counts/4;
	//	//	bool interior=(cell[1]>offset[1]&&cell[1]<lattice.cell_counts[1]-offset[1])&&
	//	//		(d<3||(cell[2]>offset[2]&&cell[2]<lattice.cell_counts[2]-offset[2]));
	//	//	elas.material_id(cell)=interior?0:1;}break;
	//	//case 5:{elas.material_id(cell)=(rand()%2==0?0:1);}break;
	//	//case 6:{
	//	//	int c=0;for(int i=0;i<d;i++){
	//	//		if(cell[i]==0||cell[i]==lattice.cell_counts[i]-1)c++;
	//	//	}
	//	//	elas.material_id(cell)=(c>=d-1?0:-1);
	//	//}
	//	//}
	//}

	//std::string file_name="D:/Bo/Data/codes/boMITLib/data/obj_grids/bunny64-3d.txt";
	//Field<int,d> voxel;
	////if(Read_Field_From_Txt_File(file_name,voxel))std::cout<<"read file "<<file_name<<std::endl;
	////else std::cerr<<"cannot read file"<<file_name<<std::endl;
	//
	//std::ifstream file(file_name,std::ios::in);if(!file){std::cout<<"cannot open file "<<file_name<<std::endl;}
	//VectorDi voxel_cell_counts=VectorDi ::Zero();for(int i=0;i<d;i++)file>>voxel_cell_counts[i];
	//voxel.Resize(voxel_cell_counts);int n=voxel_cell_counts.prod();
	//for(int i=0;i<n;i++){file>>voxel.array[i];}

	//Lattice<d> voxel_lattice(voxel.counts,(T)1);
	//iterate_cell(iter,lattice){const VectorDi& cell=iter.Coord();
	//	elas.material_id(cell)=((voxel_lattice.Valid_Cell(cell)&&voxel(cell)!=0)?0:-1);
	//}
	//elas.use_irregular_domain=true;


	std::string file_name="D:/Bo/Data/codes/boMITLib/data/obj_grids/gripper32.txt";
	Field<int,d> voxel;
	//if(Read_Field_From_Txt_File(file_name,voxel))std::cout<<"read file "<<file_name<<std::endl;
	//else std::cerr<<"cannot read file"<<file_name<<std::endl;
	
	std::ifstream file(file_name,std::ios::in);if(!file){std::cout<<"cannot open file "<<file_name<<std::endl;}
	VectorDi voxel_cell_counts=VectorDi ::Zero();for(int i=0;i<d;i++)file>>voxel_cell_counts[i];
	voxel.Resize(voxel_cell_counts);int n=voxel_cell_counts.prod();
	for(int i=0;i<n;i++){file>>voxel.array[i];}

	Lattice<d> voxel_lattice(voxel.counts,(T)1);
	iterate_cell(iter,lattice){const VectorDi& cell=iter.Coord();
		elas.material_id(cell)=((voxel_lattice.Valid_Cell(cell)&&voxel(cell)!=0)?0:-1);
	}
	elas.use_irregular_domain=true;
	elas.num_multigrid_level=3;

	elas.Initialize_DoFs();

	////boundary conditions and loads
    VectorD f=VectorD::Zero();f[0]=(T)1;/*f[1]=(T)-1;*/VectorD g=VectorD::Zero();g[1]=(T)-1;
    for(LatticeIterator<d> iter(lattice);iter.Valid();iter.Next()){
        VectorDi node=iter.Coord();int node_index=lattice.Node_Index(node);
			
			//if(node[1]==2){
			//	for(int i=0;i<d;i++)elas.psi_D.push_back(BoundaryCondition(node_index,i,0));
			//	//std::cout<<"psi_D: node: "<<node.transpose()<<std::endl;
			//}
			//VectorD f=VectorD::Unit(1)*(T).1;
			//if(node[0]==lattice.node_counts[0]/2&&node[1]==lattice.node_counts[1]/2){
			//	elas.Add_Force(node_index,f);
			//	//std::cout<<"add f: node: "<<node.transpose()<<std::endl;
			//}

			if(node[0]==2&&node[1]==2){
				for(int i=0;i<d;i++)elas.psi_D.push_back(BoundaryCondition(node_index,i,0));
				std::cout<<"psi_D: node: "<<node.transpose()<<std::endl;
			}
			
			if(node[0]==3&&node[1]==3){
				VectorD f=VectorD::Unit(1)*(T).1;
				elas.Add_Force(node_index,f);
				std::cout<<"add f: node: "<<node.transpose()<<std::endl;
			}

		//////beam under gravity
		//if(node[0]==0){for(int i=0;i<d;i++)elas.psi_D.push_back(BoundaryCondition(node_index,i,0));}
		//elas.Add_Force(node_index,g);

		//if(node[0]==lattice.node_counts[0]-1&&
		//	node[1]==0&&(d==3&&node[2]==lattice.node_counts[2]/2))elas.Add_Force(node_index,g);

		////stretching on two sides
		//if(node[0]==0)elas.Add_Force(node_index,-f);
		//else if(node[0]==lattice.node_counts[0]-1)elas.Add_Force(node_index,f);

		//VectorD offset=VectorD::Unit(0)*(T).1;
		//if(node[0]==0){for(int i=0;i<d;i++)elas.psi_D.push_back(BoundaryCondition(node_index,i,0));}
		
		//else if(node[0]==lattice.node_counts[0]-1){
		//	for(int i=0;i<d;i++)elas.psi_D.push_back(BoundaryCondition(node_index,i,offset[i]));
		//}
		//if(node[0]==0)elas.Add_Force(node_index,-f);
		//else if(node[0]==lattice.node_counts[0]-1)elas.Add_Force(node_index,f);

		//VectorD dis0=VectorD::Unit(0)*(T)-.1;
		//if(node[0]==0)for(int i=0;i<d;i++)elas.psi_D.push_back(BoundaryCondition(node_index,i,dis0[i]));
		//VectorD dis1=VectorD::Unit(0)*(T).1;
		//if(node[0]==lattice.node_counts[0]-1)for(int i=0;i<d;i++)elas.psi_D.push_back(BoundaryCondition(node_index,i,dis1[i]));

		//////Cantilever beam
		//f=VectorD::Zero();f[1]=(T)-1;
		//if(node[0]==0)for(int i=0;i<d;i++)elas.psi_D.push_back(BoundaryCondition(node_index,i,0));
		//if(node[0]==lattice.node_counts[0]-1&&node[1]==lattice.node_counts[1]-1)elas.Add_Force(node_index,f);
    }
	elas.Update_Matrices();
	elas.Solve();
    //std::cout<<"u \n"<<elas.u.transpose()<<std::endl;
    //std::cout<<"f \n"<<elas.f.transpose()<<std::endl;

	////write to files
    if(!Directory_Exists(output_dir.c_str()))Create_Directory(output_dir);
    ////write lattice
    {std::string file_name=output_dir+"/lattice";
	Write_Lattice(file_name,elas.lattice);}
    ////write displacement
    {std::string file_name=output_dir+"/displacement";
	Write_Indirect_Vector_Field(file_name,elas.u,elas.matrix_to_lattice,elas.lattice.node_counts);}
    ////write force
    {std::string file_name=output_dir+"/force";
	Write_Indirect_Vector_Field(file_name,elas.f,elas.matrix_to_lattice,elas.lattice.node_counts);}
 //   ////write strain 
 //   {std::string file_name=output_dir+(d==2?"/strains2":"/strains");
 //   Field<MatrixD,d> strains(elas.lattice.cell_counts);elas.Compute_Strain(strains);
	//Write_Cell_Field(file_name,strains.array,elas.lattice.cell_counts);}
 //   ////write stress
 //   {std::string file_name=output_dir+(d==2?"/stress2":"/stress");
 //   Field<MatrixD,d> stress(elas.lattice.cell_counts);elas.Compute_Stress(stress);
 //   Write_Cell_Field(file_name,stress.array,elas.lattice.cell_counts);}
	////write material id
	{std::string file_name=output_dir+"/material";
	Write_Cell_Field(file_name,elas.material_id.array,elas.lattice.cell_counts);}
}

void TestHexFEMPsiP(PARSE_ARGS& parse_args)
{
	std::string output_dir=parse_args.Get_String_Value("-o");output_dir+="/0";
    const int s=parse_args.Get_Integer_Value("-s");
    const int solver_type=parse_args.Get_Integer_Value("-solver");
	const int processor_type=parse_args.Get_Integer_Value("-proc");
    const bool write_residual=parse_args.Get_Option_Value("-wr");

	T length=(T)2;VectorDi cell_counts=VectorDi::Ones()*s;
	//if(d>1)cell_counts[1]=cell_counts[0]/2;if(d>2)cell_counts[2]=cell_counts[0]/2;
    Lattice<d> lattice(cell_counts,length/(T)cell_counts[0]);

	ElasticHexFEM<d> elas(lattice);
	elas.solver_type=ElasticHexFEM<d>::SolverType(solver_type);
	elas.processor_type=ElasticHexFEM<d>::ProcessorType(processor_type);
	elas.des_iter_num=parse_args.Get_Integer_Value("-n_d");
	elas.asc_iter_num=parse_args.Get_Integer_Value("-n_a");
	elas.v_cycle_iter_num=parse_args.Get_Integer_Value("-n_c");
	elas.krylov_iter_num=parse_args.Get_Integer_Value("-krylov_iter_num");
	int material_pattern=parse_args.Get_Integer_Value("-mt");
	elas.use_implicit_psi_P=true;
	elas.num_multigrid_level=1;

	////reinitialize materials
	elas.materials.clear();elas.materials.push_back(ElasticMaterial((T)1,(T).3));
	for(LatticeIterator<d> iter(lattice,LatticeIterator<d>::CELL);iter.Valid();iter.Next()){elas.material_id(iter.Coord())=0;}
	//for(LatticeIterator<d> iter(lattice,LatticeIterator<d>::CELL);iter.Valid();iter.Next()){
	//	VectorDi cell=iter.Coord();
	//	switch(material_pattern){
	//	case 0:{elas.material_id(cell)=0;}break;
	//	case 1:{elas.material_id(cell)=(cell[0]<cell_counts[0]/2?0:1);}break;
	//	case 2:{int s=16;elas.material_id(cell)=((cell[0]/s)%2==0?0:1);}break;
	//	case 3:{int s=5;elas.material_id(cell)=((cell[1]/s)%2==0?0:1);}break;
	//	case 4:{
	//		VectorDi offset=lattice.cell_counts/4;
	//		bool interior=(cell[1]>offset[1]&&cell[1]<lattice.cell_counts[1]-offset[1])&&
	//			(d<3||(cell[2]>offset[2]&&cell[2]<lattice.cell_counts[2]-offset[2]));
	//		elas.material_id(cell)=interior?0:1;}break;
	//	case 5:{elas.material_id(cell)=(rand()%2==0?0:1);}break;
	//	case 6:{
	//		int c=0;for(int i=0;i<d;i++){
	//			if(cell[i]==0||cell[i]==lattice.cell_counts[i]-1)c++;
	//		}
	//		elas.material_id(cell)=(c>=d-1?0:-1);
	//	}
	//	}
	//}

	////boundary conditions and loads
    MatrixD unit_strain=MatrixD::Zero();//unit_strain(1,1)=(T)1;
	unit_strain(0,1)=(T).5;unit_strain(1,0)=(T).5;
	elas.psi_P_strain=unit_strain;
    for(LatticeIterator<d> iter(lattice);iter.Valid();iter.Next()){
        VectorDi node=iter.Coord();int node_index=lattice.Node_Index(node);
		
		////Dirichlet boundary
		if((node[0]==0/*||node[0]==lattice.node_counts[0]-1*/)&&(node[1]==0/*||node[1]==lattice.node_counts[1]-1*/)){
			VectorD pos=lattice.Node(node)-lattice.Node(VectorDi::Zero());
			VectorD node_u=unit_strain*pos;
			std::cout<<"psi_D node "<<node.transpose()<<": "<<node_u.transpose()<<std::endl;
			for(int i=0;i<d;i++)elas.psi_D.push_back(BoundaryCondition(node_index,i,node_u[i]));
		}
		////Periodic boundary
		else if(!elas.use_implicit_psi_P&&(node[0]==0||node[1]==0)){
			int axis=node[0]==0?0:1;
			VectorDi node2=node;node2[axis]=lattice.node_counts[axis]-1;
			int node_index2=lattice.Node_Index(node2);
			VectorD offset=lattice.Node(node2)-lattice.Node(node);
			VectorD delta_u=unit_strain*offset;
			for(int i=0;i<d;i++)elas.psi_P.push_back(BoundaryCondition(node_index,i,delta_u[i],node_index2));
			std::cout<<"psi_P node "<<node.transpose()<<", "<<node2.transpose()<<": "<<delta_u.transpose()<<std::endl;
		}
    }
	std::cout<<"psi_D: "<<elas.psi_D.size()<<", psi_P: "<<elas.psi_P.size()<<std::endl;
	elas.Initialize_DoFs();
	elas.Update_Matrices();
	elas.Solve();
    //std::cout<<"u \n"<<elas.u.transpose()<<std::endl;
    //std::cout<<"f \n"<<elas.f.transpose()<<std::endl;

	////write to files
    if(!Directory_Exists(output_dir.c_str()))Create_Directory(output_dir);
    ////write lattice
    {std::string file_name=output_dir+"/lattice";
	Write_Lattice(file_name,elas.lattice);}
    ////write displacement
    {std::string file_name=output_dir+"/displacement";
	Write_Indirect_Vector_Field(file_name,elas.u,elas.matrix_to_lattice,elas.lattice.node_counts);}
    ////write force
    {std::string file_name=output_dir+"/force";
	Write_Indirect_Vector_Field(file_name,elas.f,elas.matrix_to_lattice,elas.lattice.node_counts);}
    ////write strain 
    {std::string file_name=output_dir+(d==2?"/strains2":"/strains");
    Field<MatrixD,d> strains(elas.lattice.cell_counts);elas.Compute_Strain(strains);
	Write_Cell_Field(file_name,strains.array,elas.lattice.cell_counts);}
    ////write stress
    {std::string file_name=output_dir+(d==2?"/stress2":"/stress");
    Field<MatrixD,d> stress(elas.lattice.cell_counts);elas.Compute_Stress(stress);
    Write_Cell_Field(file_name,stress.array,elas.lattice.cell_counts);}
	////write material id
	{std::string file_name=output_dir+"/material";
	Write_Cell_Field(file_name,elas.material_id.array,elas.lattice.cell_counts);}
}

void TestHomogenization(PARSE_ARGS& parse_args)
{
	std::string output_dir=parse_args.Get_String_Value("-o");output_dir+="/0";
    const int s=parse_args.Get_Integer_Value("-s");
    const int solver_type=parse_args.Get_Integer_Value("-solver");
	const int processor_type=parse_args.Get_Integer_Value("-proc");
    const bool write_residual=parse_args.Get_Option_Value("-wr");
	const int num_mg_level=parse_args.Get_Integer_Value("-mg_lv");

	T length=(T)1;VectorDi cell_counts=VectorDi::Ones()*s;
	//if(d>1)cell_counts[1]=cell_counts[0]/2;if(d>2)cell_counts[2]=cell_counts[0]/2;
    Lattice<d> lattice(cell_counts,length/(T)cell_counts[0]);
	ElasticHexFEM<d> elas(lattice);
	elas.solver_type=ElasticHexFEM<d>::SolverType(solver_type);
	elas.processor_type=ElasticHexFEM<d>::ProcessorType(processor_type);
	elas.des_iter_num=parse_args.Get_Integer_Value("-n_d");
	elas.asc_iter_num=parse_args.Get_Integer_Value("-n_a");
	elas.v_cycle_iter_num=parse_args.Get_Integer_Value("-n_c");
	elas.krylov_iter_num=parse_args.Get_Integer_Value("-krylov_iter_num");
	int material_pattern=parse_args.Get_Integer_Value("-mt");
	elas.use_implicit_psi_P=true;
	elas.num_multigrid_level=num_mg_level;

	////reinitialize materials
	elas.materials.clear();elas.materials.push_back(ElasticMaterial((T)1,(T).3));
	for(LatticeIterator<d> iter(lattice,LatticeIterator<d>::CELL);iter.Valid();iter.Next()){elas.material_id(iter.Coord())=0;}

	Homogenization<d> homo(elas);
	homo.Initialize();
	homo.Homogenize();
	
	////write to files
    if(!Directory_Exists(output_dir.c_str()))Create_Directory(output_dir);
    ////write lattice
    {std::string file_name=output_dir+"/lattice";
	Write_Lattice(file_name,elas.lattice);}
    ////write displacement
    {std::string file_name=output_dir+"/displacement";
	Write_Indirect_Vector_Field(file_name,elas.u,elas.matrix_to_lattice,elas.lattice.node_counts);}
    ////write force
    {std::string file_name=output_dir+"/force";
	Write_Indirect_Vector_Field(file_name,elas.f,elas.matrix_to_lattice,elas.lattice.node_counts);}
 //   ////write strain 
 //   {std::string file_name=output_dir+(d==2?"/strains2":"/strains");
 //   Field<MatrixD,d> strains(elas.lattice.cell_counts);elas.Compute_Strain(strains);
	//Write_Cell_Field(file_name,strains.array,elas.lattice.cell_counts);}
 //   ////write stress
 //   {std::string file_name=output_dir+(d==2?"/stress2":"/stress");
 //   Field<MatrixD,d> stress(elas.lattice.cell_counts);elas.Compute_Stress(stress);
 //   Write_Cell_Field(file_name,stress.array,elas.lattice.cell_counts);}
	////write material id
	{std::string file_name=output_dir+"/material";
	Write_Cell_Field(file_name,elas.material_id.array,elas.lattice.cell_counts);}
}

void TestFlexureFEM(PARSE_ARGS& parse_args)
{
	//std::string output_dir=parse_args.Get_String_Value("-o");output_dir+="/0";
 //   const int s=parse_args.Get_Integer_Value("-s");
 //   const int solver_type=parse_args.Get_Integer_Value("-solver");
	//const int processor_type=parse_args.Get_Integer_Value("-proc");
 //   const bool write_residual=parse_args.Get_Option_Value("-wr");

	//T diameter=(T).2;T mirror_thickness=(T).03;T flexure_thickness=(T).02;T flexure_height=(T).1;T box_size=(diameter+flexure_thickness*(T)2);
	//VectorDi cell_counts=VectorDi::Ones()*s;cell_counts[1]=(int)((flexure_height/diameter)*(T)cell_counts[0]);
 //   Lattice<d> lattice(cell_counts,box_size/(T)cell_counts[0]);
	//ElasticHexFEM<d> elas(lattice);
	//elas.solver_type=ElasticHexFEM<d>::SolverType(solver_type);
	//elas.processor_type=ElasticHexFEM<d>::ProcessorType(processor_type);
	//elas.des_iter_num=parse_args.Get_Integer_Value("-n_d");
	//elas.asc_iter_num=parse_args.Get_Integer_Value("-n_a");
	//elas.v_cycle_iter_num=parse_args.Get_Integer_Value("-n_c");
	//elas.krylov_iter_num=parse_args.Get_Integer_Value("-krylov_iter_num");
	//
	//VectorD center=(T).5*(lattice.domain_min+lattice.domain_max);T R=diameter*(T).5;T r=diameter*(T).25;

	//////reinitialize materials
	//elas.materials.clear();
	//elas.materials.push_back(ElasticMaterial((T)9.1e10,(T).24));	////for mirror
	//elas.materials.push_back(ElasticMaterial((T)9.1e9,(T).24));		////for flexure
	//Array<Box<d> > flexure_boxes;
	//Vector3 box_center(box_size*(T).5,(T)0,(T)0);
	//flexure_boxes.push_back(Box<3>(Vector3(box_size*(T).5-flexure_thickness*(T).5,(T)0,(T)0),
	//	Vector3(box_size*(T).5+flexure_thickness*(T).5,flexure_height,flexure_thickness)));
	//flexure_boxes.push_back(Box<3>(Vector3((T)0,(T)0,box_size*(T).5-flexure_thickness*(T).5),
	//	Vector3(flexure_thickness,flexure_height,box_size*(T).5+flexure_thickness*(T).5)));
	//flexure_boxes.push_back(Box<3>(Vector3(box_size-flexure_thickness,(T)0,box_size*(T).5-flexure_thickness*(T).5),
	//	Vector3(box_size,flexure_height,box_size*(T).5+flexure_thickness*(T).5)));

	//for(LatticeIterator<d> iter(lattice,LatticeIterator<d>::CELL);iter.Valid();iter.Next()){
	//	VectorDi cell=iter.Coord();VectorD pos=lattice.Center(cell);
	//	T r0=sqrt(pow(pos[0]-center[0],2)+pow(pos[2]-center[2],2));
	//	if(r0>=r&&r0<=R&&pos[1]>flexure_height-mirror_thickness){elas.material_id(cell)=0;}
	//	else{
	//		bool is_flexture=false;
	//		for(int i=0;i<(int)flexure_boxes.size();i++){
	//			if(flexure_boxes[i].Inside(pos)){is_flexture=true;break;}
	//		}
	//		elas.material_id(cell)=is_flexture?1:-1;
	//	}
	//}
	//elas.Initialize_DoFs();

	//////boundary conditions and loads
	//Field<int,d> psi_D;psi_D.Resize(lattice.node_counts);psi_D.Fill(0);
 //   VectorD f=VectorD::Zero();f[1]=(T)1;
 //   for(LatticeIterator<d> iter(lattice);iter.Valid();iter.Next()){
 //       VectorDi node=iter.Coord();int node_index=lattice.Node_Index(node);
	//	if(!lattice.Valid_Node(node,elas.material_id.array))continue;
	//	//if(node[0]==0||node[2]==0||node[2]==lattice.node_counts[2]-1){
	//	if(node[1]==0){
	//		for(int i=0;i<d;i++)elas.psi_D.push_back(BoundaryCondition(node_index,i,0));
	//		psi_D(node)=1;
	//	}
	//	else elas.Add_Force(node_index,f);
 //   }
	//elas.Update_Matrices();
	//elas.Solve();
 //   //std::cout<<"u \n"<<elas.u.transpose()<<std::endl;
 //   //std::cout<<"f \n"<<elas.f.transpose()<<std::endl;

	//////write to files
 //   if(!Directory_Exists(output_dir.c_str()))Create_Directory(output_dir);
 //   ////write lattice
 //   {std::string file_name=output_dir+"/lattice";
	//Write_Lattice(file_name,elas.lattice);}
 //   ////write displacement
 //   {std::string file_name=output_dir+"/displacement";
	//Write_Indirect_Vector_Field(file_name,elas.u,elas.matrix_to_lattice,elas.lattice.node_counts);}
 //   ////write force
 //   {std::string file_name=output_dir+"/force";
	//Write_Indirect_Vector_Field(file_name,elas.f,elas.matrix_to_lattice,elas.lattice.node_counts);}
 //   ////write strain 
 //   {std::string file_name=output_dir+(d==2?"/strains2":"/strains");
 //   Field<MatrixD,d> strains(elas.lattice.cell_counts);elas.Compute_Strain(strains);
	//Write_Cell_Field(file_name,strains.array,elas.lattice.cell_counts);}
 //   ////write stress
 //   {std::string file_name=output_dir+(d==2?"/stress2":"/stress");
 //   Field<MatrixD,d> stress(elas.lattice.cell_counts);elas.Compute_Stress(stress);
 //   Write_Cell_Field(file_name,stress.array,elas.lattice.cell_counts);}
	//////write material id
	//{std::string file_name=output_dir+"/material";
	//Write_Cell_Field(file_name,elas.material_id.array,elas.lattice.cell_counts);}
	//////write psi_D
	//{std::string file_name=output_dir+"/psi_D";
	//Write_Cell_Field(file_name,psi_D.array,lattice.node_counts);}
}

void TestTetFEM(PARSE_ARGS& parse_args)
{
	std::string output_dir=parse_args.Get_String_Value("-o");output_dir+="/0";
    const int s=parse_args.Get_Integer_Value("-s");
    const int solver_type=parse_args.Get_Integer_Value("-solver");
	const int processor_type=parse_args.Get_Integer_Value("-proc");

	T length=(T)1;VectorDi cell_counts=VectorDi::Ones();cell_counts[0]=s;
	if(d>1)cell_counts[1]=cell_counts[0]/2;if(d>2)cell_counts[2]=cell_counts[0]/2;
	VolumetricMesh<d> mesh;Initialize_Latticed_Mesh(cell_counts,(T)1/cell_counts[0],mesh);
	ElasticTetFEM<d> elas(mesh);
	std::cout<<"mesh "<<mesh.vertices.size()<<", "<<mesh.elements.size()<<std::endl;
	elas.solver_type=ElasticFEM<d>::SolverType(solver_type);
	elas.processor_type=ElasticFEM<d>::ProcessorType(processor_type);
	elas.des_iter_num=parse_args.Get_Integer_Value("-n_d");
	elas.asc_iter_num=parse_args.Get_Integer_Value("-n_a");
	elas.v_cycle_iter_num=parse_args.Get_Integer_Value("-n_c");
	elas.krylov_iter_num=parse_args.Get_Integer_Value("-krylov_iter_num");

	int axis=0;VectorD f=VectorD::Zero();f[axis]=(T)20;VectorD g=VectorD::Zero();g[1]=(T)-1;
	for(int i=0;i<(int)mesh.vertices.size();i++){
		const VectorD& pos=mesh.vertices[i];
		if(pos[axis]==0){for(int j=0;j<d;j++)elas.psi_D.push_back(BoundaryCondition(i,j,0));}
		//else if(pos[axis]==1)elas.Add_Force(i,f);
		////body force
		elas.Add_Force(i,g);
	}
	elas.Update_Matrices();
	elas.Solve();
    //std::cout<<"u \n"<<elas.u<<std::endl;
    //std::cout<<"f \n"<<elas.f<<std::endl;

    if(!Directory_Exists(output_dir.c_str()))Create_Directory(output_dir);
    {std::string file_name=output_dir+"/mesh";
	VolumetricMesh<3> mesh3;Convert_Mesh(mesh,mesh3);Write_Binary_To_File(file_name,mesh3);
    std::cout<<"write to file "<<file_name<<std::endl;}
    ////write displacement
    {std::string file_name=output_dir+"/displacement";
    Field<Vector3,1> displacements(Vector1i::Ones()*((int)mesh.vertices.size()));displacements.Fill(Vector3::Zero());
    Copy(elas.u,d,displacements.array);Write_Binary_To_File(file_name,displacements);
    std::cout<<"write to file "<<file_name<<std::endl;}
    ////write force
    {std::string file_name=output_dir+"/force";
    Field<Vector3,1> forces(Vector1i::Ones()*((int)mesh.vertices.size()));forces.Fill(Vector3::Zero());
    Copy(elas.f,d,forces.array);Write_Binary_To_File(file_name,forces);
    std::cout<<"write to file "<<file_name<<std::endl;}
	////write strain
    {std::string file_name=output_dir+(d==2?"/strains2":"/strains");
    Field<MatrixD,1> strains(Vector1i::Ones()*(int)mesh.elements.size());elas.Compute_Strain(strains.array);
    Write_Binary_To_File(file_name,strains);std::cout<<"write to file "<<file_name<<std::endl;}
    ////write stress
    {std::string file_name=output_dir+(d==2?"/stress2":"/stress");
    Field<MatrixD,1> stress(Vector1i::Ones()*(int)mesh.elements.size());elas.Compute_Stress(stress.array);
    Write_Binary_To_File(file_name,stress);std::cout<<"write to file "<<file_name<<std::endl;}
}

void TestHexTetFEM(PARSE_ARGS& parse_args)
{
	//std::string output_dir=parse_args.Get_String_Value("-o");
 //   const int s=parse_args.Get_Integer_Value("-s");
 //   const bool write_residual=parse_args.Get_Option_Value("-wr");

	//T length=(T)1;VectorDi cell_counts=VectorDi::Ones();cell_counts[0]=s;
	//if(d>1)cell_counts[1]=cell_counts[0];if(d>2)cell_counts[2]=cell_counts[0];
	////VolumetricMesh<d> mesh;Initialize_Latticed_Mesh(cell_counts,(T)1/cell_counts[0],mesh);
	//Lattice<d> lattice(cell_counts,(T)1/(T)cell_counts[0]);VolumetricMesh<d> mesh;
	//HybridLatticeMesh<d> hybrid_lattice_mesh(lattice,mesh);

	//int c0=lattice.cell_counts[1]/4;int c1=lattice.cell_counts[1]*3/4;
	//for(LatticeIterator<d> iter(lattice,LatticeIterator<d>::CELL);iter.Valid();iter.Next()){
	//	const VectorDi& cell=iter.Coord();
	//	if(cell[1]>=c0&&cell[1]<c1)hybrid_lattice_mesh.Add_Cell(cell);
	//}
	//for(LatticeIterator<d> iter(lattice,LatticeIterator<d>::CELL);iter.Valid();iter.Next()){
	//	const VectorDi& cell=iter.Coord();
	//	if(cell[1]<c0||cell[1]>=c1){
	//		VectorD offset=VectorD::Zero();
	//		if(cell[1]<lattice.cell_counts[1]/4){offset=VectorD::UnitY()*lattice.dx*(T).5*(c0-1-cell[1]);}
	//		else if(cell[1]>lattice.cell_counts[1]*3/4){offset=-VectorD::UnitY()*lattice.dx*(T).5*(cell[1]-c1);}
	//		hybrid_lattice_mesh.Add_Tet_Cell(cell,offset);
	//	}
	//}

	//std::cout<<"hybrid mesh v: "<<hybrid_lattice_mesh.vertices.size()<<", e "<<hybrid_lattice_mesh.mesh.elements.size()
	//	<<", cells: "<<hybrid_lattice_mesh.cells.size()<<std::endl;

	//ElasticHexTetFEM<d> elas(hybrid_lattice_mesh);
	//int axis=0;VectorD f=VectorD::Zero();f[axis]=(T)20;VectorD g=VectorD::Zero();g[1]=(T)-1;
	//for(int i=0;i<(int)hybrid_lattice_mesh.vertices.size();i++){
	//	const VectorD& pos=hybrid_lattice_mesh.vertices[i];
	//	if(pos[axis]==0){
	//		elas.psi_D.push_back(i);
	//	}
	//	elas.Add_Force(i,g);
	//}
	//elas.Update_Matrices();
	//elas.Solve();
 //   std::cout<<"u \n"<<elas.u<<std::endl;
 //   std::cout<<"f \n"<<elas.f<<std::endl;

 //   if(!Directory_Exists(output_dir.c_str()))Create_Directory(output_dir);
 //   {std::string file_name=output_dir+"/hybrid_lattice_mesh";
	//HybridLatticeMesh<3> mesh3(*(new Lattice<3>()),*(new VolumetricMesh<3>()));Convert_Mesh(hybrid_lattice_mesh,mesh3);
 //   Write_Binary_To_File(file_name,mesh3);
	//	std::cout<<"convert mesh v: "<<mesh3.vertices.size()<<", e "<<mesh3.mesh.elements.size()
	//	<<", cells: "<<mesh3.cells.size()<<std::endl;
 //   std::cout<<"write to file "<<file_name<<std::endl;}
 //   ////write displacement
 //   {std::string file_name=output_dir+"/displacement";
 //   Field<Vector3,1> displacements(Vector1i::Ones()*((int)mesh.vertices.size()));displacements.Fill(Vector3::Zero());
 //   Copy(elas.u,d,displacements.array);
 //   Write_Binary_To_File(file_name,displacements);
 //   std::cout<<"write to file "<<file_name<<std::endl;}
 ////   ////write force
 ////   {std::string file_name=output_dir+"/force";
 ////   Field<Vector3,1> forces(Vector1i::Ones()*((int)mesh.vertices.size()));forces.Fill(Vector3::Zero());
 ////   Copy(elas.f,d,forces.array);
 ////   Write_Binary_To_File(file_name,forces);
 ////   std::cout<<"write to file "<<file_name<<std::endl;}
	////////write strain
 ////   {std::string file_name=output_dir+(d==2?"/strains2":"/strains");
 ////   Field<MatrixD,1> strains(Vector1i::Ones()*(int)mesh.elements.size());elas.Compute_Strain(strains.array);
 ////   Write_Binary_To_File(file_name,strains);
 ////   std::cout<<"write to file "<<file_name<<std::endl;}
 ////   ////write stress
 ////   {std::string file_name=output_dir+(d==2?"/stress2":"/stress");
 ////   Field<MatrixD,1> stress(Vector1i::Ones()*(int)mesh.elements.size());elas.Compute_Stress(stress.array);
 ////   Write_Binary_To_File(file_name,stress);
 ////   std::cout<<"write to file "<<file_name<<std::endl;}
}

void TestHexFEMHeat(PARSE_ARGS& parse_args)
{
	std::string output_dir=parse_args.Get_String_Value("-o");output_dir+="/0";
    const int s=parse_args.Get_Integer_Value("-s");
    const int solver_type=parse_args.Get_Integer_Value("-solver");
	const int processor_type=parse_args.Get_Integer_Value("-proc");
    const bool write_residual=parse_args.Get_Option_Value("-wr");

	T length=(T)1;VectorDi cell_counts=VectorDi::Ones()*s;
	if(d>1)cell_counts[1]=cell_counts[0]/4;if(d>2)cell_counts[2]=cell_counts[0]/2;
    Lattice<d> lattice(cell_counts,(T)1/(T)cell_counts[0]);
	ElasticHexFEM<d> elas(lattice);
	elas.solver_type=ElasticHexFEM<d>::SolverType(solver_type);
	elas.processor_type=ElasticHexFEM<d>::ProcessorType(processor_type);
	elas.des_iter_num=parse_args.Get_Integer_Value("-n_d");
	elas.asc_iter_num=parse_args.Get_Integer_Value("-n_a");
	elas.v_cycle_iter_num=parse_args.Get_Integer_Value("-n_c");
	elas.krylov_iter_num=parse_args.Get_Integer_Value("-krylov_iter_num");
	int material_pattern=parse_args.Get_Integer_Value("-mt");
	elas.use_heat=true;

	////reinitialize materials
	elas.materials.clear();elas.material_H.clear();
	int m=16;for(int i=0;i<m;i++){
		elas.materials.push_back(ElasticMaterial((T)1,(T).3));
		T coef=(T)i/(T)(m-1);coef=abs((T)2*coef-(T)1);coef=(T)2*coef-(T)1;coef-=(T).5;
		//std::cout<<"i: "<<coef<<std::endl;
		elas.material_H.push_back(VectorD::Ones()*coef);	
	}

	//for(LatticeIterator<d> iter(lattice,LatticeIterator<d>::CELL);iter.Valid();iter.Next()){
	//	VectorDi cell=iter.Coord();VectorD pos=lattice.Center(cell);
	//	T c=(pos[0]-lattice.domain_min[0])/(lattice.domain_max[0]-lattice.domain_min[0]);
	//	int n=(int)((T)(m-1)*c);elas.material_id(iter.Coord())=n;
	//	//std::cout<<"mat: "<<pos.transpose()<<": "<<n<<std::endl;
	//}

	for(LatticeIterator<d> iter(lattice,LatticeIterator<d>::CELL);iter.Valid();iter.Next()){
		elas.material_id(iter.Coord())=0;
	}

	elas.Initialize_DoFs();

	////boundary conditions and loads
    VectorD f=VectorD::Zero();f[0]=(T)1;/*f[1]=(T)-1;*/VectorD g=VectorD::Zero();g[1]=(T)-1;
    for(LatticeIterator<d> iter(lattice);iter.Valid();iter.Next()){
        VectorDi node=iter.Coord();int node_index=lattice.Node_Index(node);
		VectorD offset=VectorD::Unit(0)*(T).1;
		if(node[0]==0){for(int i=0;i<d;i++)elas.psi_D.push_back(BoundaryCondition(node_index,i,0));}
		//else if(node[0]==lattice.node_counts[0]-1)elas.Add_Force(node_index,f);
    }
	elas.Update_Matrices();
	elas.Solve();
    //std::cout<<"u \n"<<elas.u.transpose()<<std::endl;
    //std::cout<<"f \n"<<elas.f.transpose()<<std::endl;

	////write to files
    if(!Directory_Exists(output_dir.c_str()))Create_Directory(output_dir);
    ////write lattice
    {std::string file_name=output_dir+"/lattice";
	Write_Lattice(file_name,elas.lattice);}
    ////write displacement
    {std::string file_name=output_dir+"/displacement";
	Write_Indirect_Vector_Field(file_name,elas.u,elas.matrix_to_lattice,elas.lattice.node_counts);}
    ////write force
    {std::string file_name=output_dir+"/force";
	Write_Indirect_Vector_Field(file_name,elas.f,elas.matrix_to_lattice,elas.lattice.node_counts);}
    ////write strain 
    {std::string file_name=output_dir+(d==2?"/strains2":"/strains");
    Field<MatrixD,d> strains(elas.lattice.cell_counts);elas.Compute_Strain(strains);
	Write_Cell_Field(file_name,strains.array,elas.lattice.cell_counts);}
    ////write stress
    {std::string file_name=output_dir+(d==2?"/stress2":"/stress");
    Field<MatrixD,d> stress(elas.lattice.cell_counts);elas.Compute_Stress(stress);
    Write_Cell_Field(file_name,stress.array,elas.lattice.cell_counts);}
	////write material id
	{std::string file_name=output_dir+"/material";
	Write_Cell_Field(file_name,elas.material_id.array,elas.lattice.cell_counts);}
	////write heat
	if(elas.use_heat){std::string file_name=output_dir+"/heat";
	Write_Cell_Field(file_name,elas.heat.array,elas.lattice.cell_counts);}
}

int main(int argc,char* argv[])
{
	std::cout<<"start testing elasticity fem"<<std::endl;

    ////parse arguments
    PARSE_ARGS parse_args;
    parse_args.Add_String_Argument("-o","output","output path");
    parse_args.Add_Integer_Argument("-s",8,"resolution");
    parse_args.Add_Integer_Argument("-solver",1,"solver type");
    parse_args.Add_Option_Argument("-wr","write residual");
	parse_args.Add_Integer_Argument("-n_d",1,"n descending");
	parse_args.Add_Integer_Argument("-n_a",1,"n ascending");
	parse_args.Add_Integer_Argument("-n_c",1,"n cycle");
	parse_args.Add_Integer_Argument("-krylov_iter_num",1000,"max iter num");
	parse_args.Add_Integer_Argument("-proc",1,"processor type");
	parse_args.Add_Integer_Argument("-test",1,"test");
	parse_args.Add_Integer_Argument("-mt",0,"material pattern");
	parse_args.Add_Integer_Argument("-mg_lv",4,"multigrid levels");
    parse_args.Parse(argc,argv);

	int test=parse_args.Get_Integer_Value("-test");
	switch (test){
	case 1:TestHexFEM(parse_args);break;
	case 2:TestHexFEMPsiP(parse_args);break;
	case 3:TestHomogenization(parse_args);break;
	case 4:TestTetFEM(parse_args);break;
	case 5:TestFlexureFEM(parse_args);break;
	case 6:TestHexFEMHeat(parse_args);break;
  case 7:TestMicrostructureOpt(parse_args); break;}
	//TestHexTetFEM(parse_args);
}

#endif