#ifndef TEMPLATE_MATRIX_GMP
#define TEMPLATE_MATRIX_GMP 1



namespace math{



  template< class T>
    class Vector{

  public:
    Vector();
    virtual unsigned int size() const = 0;

    virtual const T & operator[](unsigned int i) const=0;
    virtual  T & operator[](unsigned int i) =0;

    virtual ~Vector(){};
  };

  template< class T>
  class Matrix : public Vector<T> {
 
  private:
    Matrix * transposate;

  public: 
    Matrix();

    virtual ~Matrix() {
      if (transposate) delete transposate;
    }
 

    virtual unsigned int nRows() const = 0;
    
    virtual unsigned int nCols() const;


    virtual  T & get(unsigned int row, unsigned int col) = 0;

    virtual  void set(unsigned int row, unsigned int col, const  T& val) =0;

    /** return true if matrix is supposed to be a vector and not a matrix n x 1 */
    virtual bool isVector() const = 0;

    virtual void clear() = 0;
    /**
     * \brief assign a value at indice i
     */
    //  virtual void set(unsigned int i,const T & val) = 0;
    
    /**
     * \brief substract lambda * line j to line i
     */
    void subLine(unsigned int i,unsigned int j,const T & lambda){
      for(unsigned int k=0; k < nCols() ; ++k)
	set(i , k ,  get(i, k) - ( get(j , k) * lambda ) );
    }

    /**
     * \brief multiply line i by lambda
     */
    void mulLine(unsigned int i,const T & lambda){
      for(unsigned int k=0; k < nCols(); ++k)
	set(i , k,  get(i,k)  * lambda  );
    }
    
    Matrix * transpose();
  };
  
  template< class T>
    class Transpose: public Matrix<T> {
    private :
      Matrix<T> & source;
         
    public:
      Transpose(Matrix<T> & source_p)
	: source(source_p) {
      }

      unsigned int size() const{
	return source.size();
      }

    bool isVector() const {
      return false;
    }
    
      unsigned int nRows() const {
	return source.nCols();
      }
 
      T & get(unsigned int i, unsigned int j) {
	return source.get(j,i);
      }

      void set(unsigned int i,unsigned int j, const T & val){
	source.set(j,i,val);
      }

  };

  // template code.
  template<class T>
  Vector<T>::Vector() {
  }

  template<class T>
    Matrix<T>::Matrix() : 
      Vector<T>(),
      transposate(NULL)
  {
  }

  template<class T>
    Matrix<T> * Matrix<T>::transpose(){
    if (transposate == 0){
      transposate =new Transpose<T>(*this);
    }
    return transposate;
  }
  
  template<class T>
    unsigned int Matrix<T>::nCols() const{
    if(nRows() <=0) return this->size();
    return this->size() / nRows();
    }

  

}

#endif
