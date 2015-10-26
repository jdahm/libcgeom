
#include "Vtk.hh"
#include "Delaunay.hh"
#include <sstream>
#include <fstream>
#include <cassert>
#include <iomanip>

void Delaunay::WriteVtuFiles(const std::string& prefix)
// Outputs the delaunay subdivision to pvtuFileName
// Format:
// nVerts nEdges
// for i in [1,nVerts]:
//   Vert_i->x Vert_i->y
// for i in [1,nEdges]:
//   Edge_i->Org Edge_i->Dest
{

  // write the .pvtu pvtuFile
	if(!rank()){
		std::string pvtuFileName = prefix;
		pvtuFileName += ".pvtu";
		std::ofstream pvtuFile(pvtuFileName,std::ios::out);
		assert(pvtuFile.is_open());
		pvtuFile << "<VTKFile type=\"PUnstructuredGrid\">\n";
		pvtuFile << "<PUnstructuredGrid GhostLevel=\"0\">\n";
		pvtuFile << "<PPoints>\n";
		pvtuFile << "<PDataArray type=\"Float64\" Name=\"coordinates\" NumberOfComponents=\"3\" format=\"ascii\"/>\n";
		pvtuFile << "</PPoints>\n";
		pvtuFile << "<PPointData>\n";
		pvtuFile << "</PPointData>\n";
		pvtuFile << "<PCellData>\n";
		pvtuFile << "<PDataArray type=\"Int32\" Name=\"rank\" NumberOfComponents=\"1\" format=\"ascii\"/>\n";
		pvtuFile << "</PCellData>\n";
		for (unsigned i = 0; i < numProcs(); ++i)
		  pvtuFile << "<Piece Source=\"" << prefix << i << ".vtu\"/>\n";
		pvtuFile << "</PUnstructuredGrid>\n";
		pvtuFile << "</VTKFile>\n";
	  pvtuFile.close();
	}
  std::stringstream ss;
  ss << prefix << rank() << ".vtu";
  std::string vtuFileName = ss.str();

  std::stringstream buf;
  buf << "<VTKFile type=\"UnstructuredGrid\">\n";
  buf << "<UnstructuredGrid>\n";
  buf << "<Piece NumberOfPoints=\"" << pointList.size() 
      << "\" NumberOfCells=\"" << qeList.size() <<"\">\n";
  buf << "<Points>\n";
  buf << "<DataArray type=\"Float64\" Name=\"coordinates\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (Point2d* p : pointList)
    buf << std::setprecision(16) << p->x << " "
        << std::setprecision(16) << p->y << " "
        << 0.0 << std::endl;
  buf << "</DataArray>\n";
  buf << "</Points>\n";
  buf << "<Cells>\n";
  buf << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  for (QuadEdge* qe : qeList){
    Edge* e = reinterpret_cast<Edge*>(qe);
    PointList::size_type oi = LocatePointIndex(e->Org()), di = LocatePointIndex(e->Dest());
    buf << oi << " " << di << std::endl;
  }
  buf << "</DataArray>\n";
	buf << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
	for (unsigned i = 0; i < qeList.size(); ++i)
		buf << 2*(i+1) << "\n";
  buf << "</DataArray>\n";
  buf << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
	for (unsigned i = 0; i < qeList.size(); ++i)
		buf << "3\n";
  buf << "</DataArray>\n";
  buf << "</Cells>\n";
  buf << "<PointData>\n";
	buf << "</PointData>\n";
	buf << "<CellData>\n";
	buf << "<DataArray type=\"Int32\" Name=\"rank\" NumberOfComponents=\"1\" format=\"ascii\">\n";
	for (unsigned i = 0; i < qeList.size(); ++i)
		buf << rank() << "\n";
  buf << "</DataArray>\n";
	buf << "</CellData>\n";
  buf << "</Piece>\n";
	buf << "</UnstructuredGrid>\n";
	buf << "</VTKFile>\n";
  std::ofstream vtuFile(vtuFileName,std::ios::out);
  assert(vtuFile.is_open());
  vtuFile << buf.rdbuf();
  vtuFile.close();
}