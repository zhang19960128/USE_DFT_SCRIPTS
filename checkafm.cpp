#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <sstream>
#include <algorithm>
int checkneighbor(std::vector<std::vector<double> > &posit_spin,std::vector<std::vector<double> > &posit_spin_original,int site){
    int period=4;
    double distance=0.0;
    size_t element=posit_spin.size();
    int yes=1;
    double temp=0.0;
    for(int i=0;i<element;i++){
    if(i==site){
        continue;
    }
    distance=0.0;
    for(size_t j=0;j<3;j++){
        temp=posit_spin[i][j]-posit_spin[site][j];
        temp=period*(temp/period-round(temp/period));
       distance=distance+temp*temp;
       }
    distance=sqrt(distance);
    //std::cout<<"site: "<<site<<" element: "<<i<<" distance: "<<distance<<std::endl;
    if(std::fabs(distance-2.0)<0.0001){
        yes=yes*(posit_spin[site][3]*posit_spin[i][3]<-0.0000001);
    }
    else{
    }
    if(yes==0){
    std::cout<<"not AFM found:"<<std::endl;
    for(size_t k=0;k<4;k++){
        std::cout<<posit_spin_original[i][k]<<"\t";
    }
    for(size_t k=0;k<4;k++){
        std::cout<<posit_spin_original[i][k]<<"\t";
    }
    break;
    }
    }
    return yes;
}
int main(int argc,char* argv[]){
    std::fstream fs;
    fs.open(argv[1],std::fstream::in);
    std::string tempstring;
    std::map<int,double> magmap;
    std::vector<std::string> vec_species;
    int indexnumber,index;
    std::string indexnumber_string;
    std::stringstream temp_stream;
    std::string species;
    double magnetization;
    int natom;
    double period[3][3]={{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
    std::vector<std::vector<double> > posit_spin(0,std::vector<double>(4,0.0));
    std::vector<double> work_vec(4,0.0);
    int cartesian_coord;
    while(getline(fs,tempstring)){
    if(tempstring.find("starting_magnetization")!=std::string::npos){
        indexnumber=tempstring.find("(");
        temp_stream.clear();
        temp_stream.str(tempstring.substr(indexnumber+1,1));
        temp_stream>>index;
        temp_stream.clear();
        indexnumber=tempstring.find("=");
        temp_stream.str(tempstring.substr(indexnumber+1,100));
        temp_stream>>magnetization;
        magmap.insert(std::pair<int,double>(index,magnetization));
    }
    if(tempstring.find("upf")!=std::string::npos || tempstring.find("UPF")!=std::string::npos){
    temp_stream.clear();
    temp_stream.str(tempstring);
    temp_stream>> species;
    vec_species.push_back(species);
    }
    if(tempstring.find("nat")!=std::string::npos){
    temp_stream.clear();
    indexnumber=tempstring.find("=");
    temp_stream.clear();
    temp_stream.str(tempstring.substr(indexnumber+1,100));
    temp_stream>>natom;
    temp_stream.clear();
    }
    if(tempstring.find("CELL_PARAMETERS")!=std::string::npos){
    for(size_t i=0;i<3;i++){
    getline(fs,tempstring);
    temp_stream.clear();
    temp_stream.str(tempstring);
    for(size_t j=0;j<3;j++){
    temp_stream>>period[i][j];
    }
    temp_stream.clear();
    }
    }
    if(tempstring.find("ATOMIC_POSITIONS")!=std::string::npos){
    if(tempstring.find("angstrom")!=std::string::npos){
        cartesian_coord=1;
    }
    else{
        cartesian_coord=0;
    }
    for(size_t i=0;i<natom;i++){
        getline(fs,tempstring);
        temp_stream.clear();
        temp_stream.str(tempstring);
        temp_stream>>species;
        index=0;
        for(size_t i=0;i<vec_species.size();i++){
            if(vec_species[i].find(species)!=std::string::npos){
                index=i;
                break;
            }
        }
        magnetization=magmap[index+1];
        work_vec[3]=magnetization;
        for(size_t i=0;i<3;i++){
            temp_stream>>work_vec[i];
        }
        if(std::fabs(magnetization)>0.000001){
            posit_spin.push_back(work_vec);
        }
    }
    }
    }
    fs.close();
    size_t vec_size=posit_spin.size();
    std::vector<std::vector<double> > posit_spin_original(posit_spin);
    if(cartesian_coord){
        for(size_t i=0;i<vec_size;i++){
            for(size_t j=0;j<3;j++){
                posit_spin[i][j]=round(posit_spin[i][j]/period[j][j]*4);
            }
        }
    }
    else{
        for(size_t i=0;i<vec_size;i++){
            for(size_t j=0;j<3;j++){
                posit_spin[i][j]=round(posit_spin[i][j]*4);
            }
        }
    }
    for(int i=0;i<vec_size;i++){
    std::cout<<checkneighbor(posit_spin,posit_spin_original,i)<<std::endl;
    }
}
