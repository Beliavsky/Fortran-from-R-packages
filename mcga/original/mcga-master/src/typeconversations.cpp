#include <Rcpp.h>
#include <limits>

using namespace Rcpp;

// [[Rcpp::export]]
double MaxDouble(){
  return (std::numeric_limits<double>::max());
}

// [[Rcpp::export]]
int SizeOfDouble(){
  return(sizeof(double));
}

// [[Rcpp::export]]
int SizeOfInt(){
  return(sizeof(int));
}

// [[Rcpp::export]]
int SizeOfLong(){
  return(sizeof(long));
}

// [[Rcpp::export]]
IntegerVector DoubleToBytes(double x) {
  unsigned char *c = (unsigned char*)&x;
  IntegerVector result(sizeof(double));
  for (unsigned int i=0; i<sizeof(double); i++){
    result[i] = *c;
    c++;
  }
  return(result);
}


// [[Rcpp::export]]
IntegerVector DoubleVectorToBytes(DoubleVector d){
  int a = 0;
  IntegerVector byte_result(d.length() * sizeof(double));
  for (int i=0;i < d.length(); i++){
    double current = d[i];
    IntegerVector v = DoubleToBytes(current);
    for (unsigned int j=0; j < sizeof(double); j++ ){
      byte_result[a] = v[j];
      a++;
    }
  }
  return(byte_result);
}

// [[Rcpp::export]]
double BytesToDouble(IntegerVector x){
  double d;
  unsigned char *c = (unsigned char*)&d;
  for (unsigned int i=0; i<sizeof(double); i++){
    *c = x[i];
    c++;
  }
  return(d);
}

// [[Rcpp::export]]
DoubleVector ByteVectorToDoubles(IntegerVector b){
  unsigned int numDoubles = (unsigned int) (b.length() / sizeof(double));
  DoubleVector doubles(numDoubles);
  unsigned int a = 0;
  for (unsigned int i = 0; i < numDoubles; i++){
    IntegerVector currentBytes(sizeof(double));
    for (unsigned int j=0; j < sizeof(double); j++){
      currentBytes[j] = b[a];
      a++;
    }
    doubles[i] = BytesToDouble(currentBytes);
  }
  return(doubles);
}


// [[Rcpp::export]]
List OnePointCrossOver(IntegerVector bytes1, IntegerVector bytes2, unsigned int cutpoint){
  unsigned int len = bytes1.length();
  IntegerVector bytes3(len);
  IntegerVector bytes4(len);
  for (unsigned int i=0; i< len; i++){
    if(i < cutpoint){
      bytes3[i] = bytes1[i];
      bytes4[i] = bytes2[i];
    }else{
      bytes3[i] = bytes2[i];
      bytes4[i] = bytes1[i];
    }
  }
  List out(2);
  out[0] = bytes3;
  out[1] = bytes4;
  return(out);
}


// [[Rcpp::export]]
List OnePointCrossOverOnDoublesUsingBytes(DoubleVector d1, DoubleVector d2, unsigned int cutpoint){
  IntegerVector bytes1 = DoubleVectorToBytes(d1);
  IntegerVector bytes2 = DoubleVectorToBytes(d2);
  unsigned int len = bytes1.length();
  IntegerVector bytes3(len);
  IntegerVector bytes4(len);
  for (unsigned int i=0; i< len; i++){
    if(i < cutpoint){
      bytes3[i] = bytes1[i];
      bytes4[i] = bytes2[i];
    }else{
      bytes3[i] = bytes2[i];
      bytes4[i] = bytes1[i];
    }
  }
  List out(2);
  out[0] = ByteVectorToDoubles(bytes3);
  out[1] = ByteVectorToDoubles(bytes4);
  return(out);
}


// [[Rcpp::export]]
List TwoPointCrossOver(IntegerVector bytes1, IntegerVector bytes2, unsigned int cutpoint1, unsigned int cutpoint2){
  unsigned int len = bytes1.length();
  IntegerVector bytes3(len);
  IntegerVector bytes4(len);
  for (unsigned int i=0; i< len; i++){
    if( (i < cutpoint1) || (i > cutpoint2) ){
      bytes3[i] = bytes1[i];
      bytes4[i] = bytes2[i];
    }else{
      bytes3[i] = bytes2[i];
      bytes4[i] = bytes1[i];
    }
  }
  List out(2);
  out[0] = bytes3;
  out[1] = bytes4;
  return(out);
}


// [[Rcpp::export]]
List TwoPointCrossOverOnDoublesUsingBytes(DoubleVector d1, DoubleVector d2, unsigned int cutpoint1, unsigned int cutpoint2){
  IntegerVector bytes1 = DoubleVectorToBytes(d1);
  IntegerVector bytes2 = DoubleVectorToBytes(d2);
  unsigned int len = bytes1.length();
  IntegerVector bytes3(len);
  IntegerVector bytes4(len);
  DoubleVector off1;
  DoubleVector off2;
  for (unsigned int i=0; i< len; i++){
    if( (i < cutpoint1) || (i > cutpoint2) ){
      bytes3[i] = bytes1[i];
      bytes4[i] = bytes2[i];
    }else{
      bytes3[i] = bytes2[i];
      bytes4[i] = bytes1[i];
    }
  }
  off1 = ByteVectorToDoubles(bytes3);
  off2 = ByteVectorToDoubles(bytes4);
  List out(2);
  out[0] = off1;
  out[1] = off2;
  return(out);
}


// [[Rcpp::export]]
List UniformCrossOver(IntegerVector bytes1, IntegerVector bytes2){
  unsigned int len = bytes1.length();
  IntegerVector bytes3(len);
  IntegerVector bytes4(len);
  DoubleVector u = runif(len, 0.0, 1.0);
  for (unsigned int i=0; i< len; i++){
    if(u[i] < 0.5){
      bytes3[i] = bytes1[i];
      bytes4[i] = bytes2[i];
    }else{
      bytes3[i] = bytes2[i];
      bytes4[i] = bytes1[i];
    }
  }
  List out(2);
  out[0] = bytes3;
  out[1] = bytes4;
  return(out);
}


// [[Rcpp::export]]
List UniformCrossOverOnDoublesUsingBytes(DoubleVector d1, DoubleVector d2){
  unsigned int len;
  IntegerVector bytes1 = DoubleVectorToBytes(d1);
  IntegerVector bytes2 = DoubleVectorToBytes(d2);
  len = bytes1.length();
  IntegerVector bytes3(len);
  IntegerVector bytes4(len);
  DoubleVector off1(len);
  DoubleVector off2(len);
  DoubleVector u = runif(len, 0.0, 1.0);
  
  for (unsigned int i=0; i< len; i++){
    if(u[i] < 0.5){
      bytes3[i] = bytes1[i];
      bytes4[i] = bytes2[i];
    }else{
      bytes3[i] = bytes2[i];
      bytes4[i] = bytes1[i];
    }
  }
  off1 = ByteVectorToDoubles(bytes3);
  off2 = ByteVectorToDoubles(bytes4);
  
  List out(2);
  out[0] = off1;
  out[1] = off2;
  return(out);
}




// [[Rcpp::export]]
IntegerVector ByteCodeMutation(IntegerVector bytes1, const double pmutation){
  unsigned int len = bytes1.length();
  DoubleVector u = runif(len, 0.0, 1.0);
  DoubleVector u2 = runif(len, 0.0, 1.0);
  IntegerVector newbytes(len);
  for (unsigned int i=0; i<len; i++){
    newbytes[i] = bytes1[i];
    if(u[i] < pmutation){
       if(u2[i] < 0.5){
         newbytes[i]++;
       }else{
         newbytes[i]--;
       }
       if(newbytes[i]<0){
         newbytes[i] = 255;
       }
       if(newbytes[i]>255){
         newbytes[i] = 0;
       }
    }
  }
  return(newbytes);
}


// [[Rcpp::export]]
DoubleVector ByteCodeMutationUsingDoubles(DoubleVector d, const double pmutation){
  IntegerVector bytes1 = DoubleVectorToBytes(d);
  unsigned int len = bytes1.length();
  DoubleVector u = runif(len, 0.0, 1.0);
  DoubleVector u2 = runif(len, 0.0, 1.0);
  IntegerVector newbytes(len);
  DoubleVector result;
  for (unsigned int i=0; i < len; i++){
    newbytes[i] = bytes1[i];
    if(u[i] < pmutation){
      if(u2[i] < 0.5){
        newbytes[i]++;
      }else{
        newbytes[i]--;
      }
      if(newbytes[i]<0){
        newbytes[i] = 255;
      }
      if(newbytes[i]>255){
        newbytes[i] = 0;
      }
    }
  }
  result = ByteVectorToDoubles(newbytes);
  return(result);
}

// [[Rcpp::export]]
DoubleVector ByteCodeMutationUsingDoublesRandom(DoubleVector d, const double pmutation){
  IntegerVector bytes1 = DoubleVectorToBytes(d);
  unsigned int len = bytes1.length();
  DoubleVector u = runif(len, 0.0, 1.0);
  IntegerVector newbytes(len);
  DoubleVector result;
  for (unsigned int i=0; i < len; i++){
    newbytes[i] = bytes1[i];
    if(u[i] < pmutation){
		newbytes[i] = round(runif(1, 0.0, 255.0)[0]);
    }
  }
  result = ByteVectorToDoubles(newbytes);
  return(result);
}


// [[Rcpp::export]]
void EnsureBounds (DoubleVector doubles, DoubleVector mins, DoubleVector maxs){
  for (int i = 0; i < doubles.length(); i++){
    if( (doubles[i] < mins[i]) || (doubles[i] > maxs[i])){
      doubles[i] = runif(1, mins[i], maxs[i])[0];
    }
  }
}



