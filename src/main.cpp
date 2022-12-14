#include <iostream>
#include "../common/mesh_io.h"
#include "../alg/deformation.h"
#include "../alg/remeshing.h"
#include "../interface/window.h"

int main(int argc, char* argv[])
{
	if (argc < 2) return -1;

	gui::Window window("window");
	Eigen::Matrix3Xd V;
	Eigen::Matrix3Xi F;
	io::read_triangle_obj(argv[1], V, F);

	Eigen::Vector3d dir;
	dir << 5, 0, 0;

	alg::Deformer er(V, F); 
	//er.mesh_deformation(dir, 22, 0.3, 0.8);
	er.mesh_deformation(dir, 8880,1, 0.8);

	er.get_result_mesh(V, F);
	io::save_triangle_obj(argv[2], V, F);

	//win.show();


	return 1;

}