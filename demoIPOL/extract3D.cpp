#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"

#include "DGtal/images/ImageSelector.h"
#include "DGtal/helpers/StdDefs.h"

#include "DGtal/io/Display3D.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/Color.h"

#include "DGtal/topology/KhalimskySpaceND.h"
#include "DGtal/topology/helpers/Surfaces.h"

#include "ImaGene/Arguments.h"




using namespace std;
using namespace DGtal;
using namespace Z3i;



static ImaGene::Arguments args;


int main( int argc, char** argv )
{


  args.addOption("-image", "-image <filename>  ", "aFile.vol ");
  args.addOption("-output", "-output <filename> the output filename with .off extension", "output.off"); 
  args.addOption("-exportSRC", "-exportSRC <filename> export the source set of voxels", "src.off"); 
  args.addOption("-threshold", "-threshold <min> <max> (default: min = 128, max 255  ", "128", "255");
  args.addOption( "-badj", "-badj <0/1>: 0 is interior bel adjacency, 1 is exterior (def. is 0).", "0" );

  if ( ( argc <= 1 ) ||  ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "extract3D: ", 
			  "Extracts all 3D connected components from a .vol 3D image and generate a resulting 3D mesh on .OFF format. \nTypical use: \n extract3D -threshold 200 -image image.pgm > imageContour.fc ",
			  "" )
	   << endl;
      return 1;
    }  
  
  string imageFileName = args.getOption("-image")->getValue(0);  
  string outputFileName = args.getOption("-output")->getValue(0);  
  string srcFileName = args.getOption("-exportSRC")->getValue(0);  
  
  int valMin = args.getOption("-threshold")->getIntValue(0);
  int valMax = args.getOption("-threshold")->getIntValue(1);
  bool badj = (args.getOption("-badj")->getIntValue(0))!=1;
  
  
  
  
  typedef ImageSelector < Domain, int>::Type Image;
  Image image =   VolReader<Image>::importVol(imageFileName);
  Z3i::DigitalSet diamond_set(image.domain());
  SetFromImage<Z3i::DigitalSet>::append<Image>(diamond_set, image, valMin, valMax);
  //A KhalimskySpace is constructed from the domain boundary points.
  Point pUpper = image.domain().upperBound();
  Point pLower = image.domain().lowerBound();
  cerr << "diamond size" << diamond_set.size(); 
  //Increase space to process also cell in the image border:
  pLower[0]-=1;pLower[1]-=1;pLower[2]-=1;
  pUpper[0]+=1;pUpper[1]+=1;pUpper[2]+=1;

  KSpace K;
  K.init(pLower, pUpper, true);
  
  SurfelAdjacency<3> sAdj(  badj );

  vector<vector<SCell> > vectConnectedSCell;
  
  
  //Here since the last argument is set to true, the resulting
  //SignedKhalimskySpaceND are signed in order to indicate the direction
  //of exterior. You can also get the SignefKhalimskySpaceND with default
  //sign:

  SetPredicate<DigitalSet> shape_set_predicate( diamond_set );
  Surfaces<KSpace>::extractAllConnectedSCell(vectConnectedSCell,K, sAdj, shape_set_predicate, false);


  
  Display3D exportSRC;
  Display3D exportSurfel;

    

  // Each connected compoments are simply displayed with a specific color.
  GradientColorMap<long> gradient(0, (const long)vectConnectedSCell.size());
  gradient.addColor(Color::Red);
  gradient.addColor(Color::Yellow);
  gradient.addColor(Color::Green);
  gradient.addColor(Color::Cyan);
  gradient.addColor(Color::Blue);
  gradient.addColor(Color::Magenta);
  gradient.addColor(Color::Red);  
 
  // Processing ViewerInt  
  for(int i=0; i< vectConnectedSCell.size();i++){
    DGtal::Color col= gradient(i);
    exportSurfel << CustomColors3D(Color(250, 0,0), Color(col.red(), 
						       col.green(),
						       col.blue()));
    for(int j=0; j< vectConnectedSCell.at(i).size();j++){
      exportSurfel << vectConnectedSCell.at(i).at(j);
    }    
  }
  
  exportSurfel << CustomColors3D(Color(250, 0,0),Color(250, 200,200, 200));
  exportSurfel << diamond_set;

  exportSRC << diamond_set;

  
  exportSurfel >> outputFileName;
  exportSRC >> srcFileName;
  
  
  
}
//                                                                           //
