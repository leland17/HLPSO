using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Diagnostics;


namespace MatrixHandle
{
    public class Matrix
    {

        //构造一维行向量
        public Matrix(int col)
        {
            m_data = new double[1, col];
            //for (int i = 0; i < col; i++)
            //{
            //    m_data[0,i] = 0; 
            //} 

        }
        //构造一般矩阵
        public Matrix(int row, int col)
        {
            m_data = new double[row, col];
            //for (int i = 0; i < row; row++)
            //{
            //    for (int j = 0; j < col; col++)
            //    {
            //        m_data[i, j] = 0;
            //    }
            //}
        }
        //复制构造函数
        public Matrix(Matrix m)
        {
            int row = m.Row;
            int col = m.Col;
            m_data = new double[row, col];

            for (int i = 0; i < row; i++)
                for (int j = 0; j < col; j++)
                    m_data[i, j] = m[i, j];

        }

       //分配方阵的大小 //对于已含有内存的矩阵，将清空数据
        public void SetSize(int row)
        {
            m_data = new double[row,row];
        }

        //分配矩阵的大小//对于已含有内存的矩阵，将清空数据
        public void SetSize(int row,int col)
        {
            m_data = new double[row,col];
        }
        
        //unit matrix:设为单位阵
        public void SetUnit()
        {
            for (int i = 0; i < m_data.GetLength(0); i++)
                for (int j = 0; j < m_data.GetLength(1); j++)
                    m_data[i, j] = ((i == j) ? 1 : 0);
        }

        //设置元素值
        public void SetMatrixValue(int i, int j,double d)
        {
               m_data[i, j] = d;
        }

        public void SetMatrixValue(int i, double d)
        {
            if (Row == 1 || Col == 1)
            {
                if (Row == 1)
                {
                    m_data[0, i] = d;
                }
                else
                {
                    m_data[i, 0] = d;
                }
            }
            else 
            {
                throw new Exception("matrix 为一阶行矩阵或一阶列矩阵");
            }
            

        }
       
        //获取元素值
        public double GetMatrixValue(int i, int j)
        {
            return m_data[i, j];
        }

        public double GetMatrixValue(int i)
        {
            if (Row == 1 || Col == 1)
            {
                if (Row == 1)
                {
                    return m_data[0, i];
                }
                else
                {
                    return m_data[i, 0];
                }

            }
            else
            {
                throw new Exception("matrix 为一阶行矩阵或一阶列矩阵");
            }

         
            
        }
        // Value extraction：返中行数
        public int Row
        {
            get
            {
                return m_data.GetLength(0);
            }
        }

        //返回列数
        public int Col
        {
            get
            {
                return m_data.GetLength(1);
            }
        }

        //重载索引 //存取数据成员
        public double this[int row, int col]
        {
            get
            {
                return m_data[row, col];
            }
            set
            {
                m_data[row, col] = value;
            }
        }

        //初等变换　对调两行：ri<-->rj
        public Matrix Exchange(int i, int j)
        {
            double temp;

            for (int k = 0; k < Col; k++)
            {
                temp = m_data[i, k];
                m_data[i, k] = m_data[j, k];
                m_data[j, k] = temp;
            }
            return this;
        }


        //初等变换　第index 行乘以mul
        Matrix Multiple(int index, double mul)
        {
            for (int j = 0; j < Col; j++)
            {
                m_data[index, j] *= mul;
            }
            return this;
        }


        //初等变换 第src行乘以mul加到第index行
        Matrix MultipleAdd(int index, int src, double mul)
        {
            for (int j = 0; j < Col; j++)
            {
                m_data[index, j] += m_data[src, j] * mul;
            }

            return this;
        }

        //transpose 转置
        public Matrix Transpose()
        {
            Matrix ret = new Matrix(Col, Row);

            for (int i = 0; i < Row; i++)
                for (int j = 0; j < Col; j++)
                {
                    ret[j,i] = m_data[i, j];
                }
            return ret;
        }

        //binary addition 矩阵加
        public static Matrix operator +(Matrix lhs, Matrix rhs)
        {
            if (lhs.Row != rhs.Row)    //异常
            {
                System.Exception e = new Exception("相加的两个矩阵的行数不等");
                throw e;
            }
            if (lhs.Col != rhs.Col)     //异常
            {
                System.Exception e = new Exception("相加的两个矩阵的列数不等");
                throw e;
            }

            int row = lhs.Row;
            int col = lhs.Col;
            Matrix ret = new Matrix(row, col);

            for (int i = 0; i < row; i++)
                for (int j = 0; j < col; j++)
                {
                    double d = lhs[i, j] + rhs[i, j];
                    ret[i, j] = d;
                }
            return ret;

        }

        //binary subtraction 矩阵减
        public static Matrix operator -(Matrix lhs, Matrix rhs)
        {
            if (lhs.Row != rhs.Row)    //异常
            {
                System.Exception e = new Exception("相减的两个矩阵的行数不等");
                throw e;
            }
            if (lhs.Col != rhs.Col)     //异常
            {
                System.Exception e = new Exception("相减的两个矩阵的列数不等");
                throw e;
            }

            int row = lhs.Row;
            int col = lhs.Col;
            Matrix ret = new Matrix(row, col);

            for (int i = 0; i < row; i++)
                for (int j = 0; j < col; j++)
                {
                    double d = lhs[i, j] - rhs[i, j];
                    ret[i, j] = d;
                }
            return ret;
        }


        //binary multiple 矩阵乘
        public static Matrix operator *(Matrix lhs, Matrix rhs)
        {
            if (lhs.Col != rhs.Row)    //异常
            {
                System.Exception e = new Exception("相乘的两个矩阵的行列数不匹配");
                throw e;
            }
            Matrix ret = new Matrix(lhs.Row, rhs.Col);
            double temp;
            for (int i = 0; i < lhs.Row; i++)
            {
                for (int j = 0; j < rhs.Col; j++)
                {
                    temp = 0;
                    for (int k = 0; k < lhs.Col; k++)
                    {
                        temp += lhs[i, k] * rhs[k, j];
                    }
                    ret[i, j] = temp;
                }
            }

            return ret;
        }


        //binary division 矩阵除
        public static Matrix operator /(Matrix lhs, Matrix rhs)
        {
            return lhs * rhs.InverseMatrix();
        }

        //unary addition单目加
        public static Matrix operator +(Matrix m)
        {
            Matrix ret = new Matrix(m);
            return ret;
        }

        //unary subtraction 单目减
        public static Matrix operator -(Matrix m)
        {
            Matrix ret = new Matrix(m);
            for (int i = 0; i < ret.Row; i++)
                for (int j = 0; j < ret.Col; j++)
                {
                    ret[i, j] = -ret[i, j];
                }

            return ret;
        }

        //number multiple 数乘
        public static Matrix operator *(double d, Matrix m)
        {
            Matrix ret = new Matrix(m);
            for (int i = 0; i < ret.Row; i++)
                for (int j = 0; j < ret.Col; j++)
                    ret[i, j] *= d;

            return ret;
        }

        //number division 数除
        public static Matrix operator /(double d, Matrix m)
        {
            return d * m.InverseMatrix();
        }

        //功能：返回列主元素的行号
        //参数：row为开始查找的行号
        //说明：在行号[row,Col)范围内查找第row列中绝对值最大的元素，返回所在行号
        int Pivot(int row)
        {
            int index = row;

            for (int i = row + 1; i < Row; i++)
            {
                if (m_data[i, row] > m_data[index, row])
                    index = i;
            }

            return index;
        }

        //求矩阵行列式，只有方阵能够求
        public double Determinant()
        {
            //二阶及以下行列式直接计算
            if (Row == 0) return 0;
            else if (Row == 1) return m_data[0,0];
            else if (Row == 2)
            {
                return m_data[0,0] * m_data[1,1] - m_data[0,1] * m_data[1,0];
            }

            //对第一行使用“加边法”递归计算行列式的值
            double dSum = 0, dSign = 1;
            for (int i = 0; i < Row; i++)
            {
                Matrix matrixTemp = new Matrix(Row - 1,Row - 1);
                
                for (int j = 0; j < matrixTemp.Row; j++)
                {
                    for (int k = 0; k < matrixTemp.Row; k++)
                    {
                        matrixTemp[j,k] = m_data[j + 1,k >= i ? k + 1 : k];
                    }
                }
                dSum += (m_data[0,i] * dSign * matrixTemp.Determinant());
                dSign = dSign * -1;
            }

            return dSum;
        }
   

        //求解伴随矩阵
        public Matrix AdjointMatrix()
        {
            //制作一个伴随矩阵大小的矩阵
            Matrix result = new Matrix(Row,Col);
            

            //生成伴随矩阵
            for (int i = 0; i < result.Row; i++)
            {
                for (int j = 0; j < result.Row; j++)
                {
                    //存储代数余子式的矩阵（行、列数都比原矩阵少1）
                    if (result.Row == 1)
                    {
                        result[0,0] = 1;
                    }
                    else 
                    {
                        Matrix temp = new Matrix(result.Row - 1, result.Row - 1);
                        //生成代数余子式
                        for (int x = 0; x < temp.Row; x++)
                        {
                            for (int y = 0; y < temp.Row; y++)
                            {
                                temp[x, y] = m_data[x < i ? x : x + 1, y < j ? y : y + 1];
                            }
                        }
                        //Console.WriteLine("代数余子式:");
                        //PrintMatrix(temp);
                        result[j, i] = ((i + j) % 2 == 0 ? 1 : -1) * temp.Determinant();
                    }
                }
            }
            //Console.WriteLine("伴随矩阵：");
            //PrintMatrix(result);
            return result;
        }
        
        //get the  colum primary element of matrix kth row,k means the kth row of the number.
        public  void selectMainElement(int k)     
        {  
          // 寻找第k列的主元素以及它所在的行号  
          double t, mainElement;            // mainElement用于保存主元素的值  
          int l;                            // 用于保存主元素所在的行号  
 
          // 从第k行到第n行寻找第k列的主元素，记下主元素mainElement和所在的行号l  
          mainElement = Math.Abs(m_data[k, k]);  // 注意别忘了取绝对值  
          l = k;  
          for(int i = k + 1; i < Row; i++)  
          {  
            if (mainElement < Math.Abs(m_data[i, k]))  
            {  
              mainElement = Math.Abs(m_data[i, k]);  
              l = i;                        // 记下主元素所在的行号  
            }  
          }  
 
          // l是主元素所在的行。将l行与k行交换，每行前面的k个元素都是0，不必交换  
          if (l != k)  
          {  
            for (int j = k; j <= Row; j++)  
            {   
              t = m_data[k, j]; m_data[k, j] = m_data[l, j]; m_data[l, j] = t;  
            }  
          }  
        }

        // get the Matrix after the Gaussian elimination;
        public Matrix Gauss()
        {
            double d;
            Matrix GaussMatrix = new Matrix(Row,Col);
           
            // 消元  
            for (int k = 0; k < Row; k++)
            {
                selectMainElement(k); // 选择主元素  
           
                // for (int j = k; j <= n; j++ ) a[k, j] = a[k, j] / a[k, k];  
                // 若将下面两个语句改为本语句，则程序会出错，因为经过第1次循环  
                // 后a[k,k]=1，a[k,k]的值发生了变化，所以在下面的语句中先用d  
                // 将a[k,k]的值保存下来  
                d = m_data[k, k];
                int m = k, n = k;
                while(-1E-9 < d && d < 1E-9)
                {
                    d = m_data[m,n];
                    n++;
                }
                for (int j = n; j <= Row; j++)
                {
                    {
                        m_data[k, j] = m_data[k, j] / d;
                    }
                } 
               
                // Guass消去法与Jordan消去法的主要区别就是在这一步，Gauss消去法是从k+1  
                // 到n循环，而Jordan消去法是从1到n循环，中间跳过第k行  
                for (int i = k + 1; i < Row; i++)
                {
                    d = m_data[i, k];  // 这里使用变量d将a[i,k]的值保存下来的原理与上面注释中说明的一样  
                    for (int j = k; j <= Row; j++)
                    {
                        m_data[i, j] = m_data[i, j] - d * m_data[k, j];
                    } 
                }
                
            }
            GaussMatrix.m_data = m_data;
            return GaussMatrix;
        }  

        //get the index of linear independent square matrix in the Guass matrix
        public int[] IndexofIndependentMatrix()
        {
            int[] Index= new int[Col];
           
            for (int i = 0; i < Row; i++)
            {
                for (int j = 0; j < Col; j++)
                {
                    if (m_data[i, j] != 0)
                    {
                        Index[j] = 1;
                        break;
                    }
                }
            }
            return Index;
        }



        //inversion 逆阵：使用矩阵的初等变换，列主元素消去法
        public Matrix InverseMatrix()
        {
            //matrix必须为非空
            if (m_data == null || Row == 0)
            {
                return new Matrix(0,0);
            }

            //matrix 必须为方阵
            int len = Row;
            for (int counter = 0; counter < Row; counter++)
            {
                if (Row != Col)
                {
                    throw new Exception("matrix 必须为方阵");
                }
            }

            //计算矩阵行列式的值
            Matrix temp = new Matrix(Row, Col);
            for (int i = 0; i < Row; i++)
            {
                for (int j = 0; j < Col; j++)
                {
                    temp.SetMatrixValue(i, j, m_data[i, j]);
                }
            }

            double dDeterminant = temp.Determinant();
            if (Math.Abs(dDeterminant) <= 1E-8)
            {
                throw new Exception("矩阵不可逆");
            }

            //制作一个伴随矩阵大小的矩阵

            Matrix result = temp.AdjointMatrix();

            //伴随矩阵的每项除以矩阵行列式的值，即为所求
            for (int i = 0; i < Row; i++)
            {
                for (int j = 0; j < Row; j++)
                {
                    result[i,j] = result[i,j] / dDeterminant;
                }
            }

            return result;
        }
       

        //determine if the matrix is square:方阵
        public bool IsSquare()
        {
            return Row == Col;
        }

        //determine if the matrix is symmetric对称阵
        public bool IsSymmetric()
        {

            if (Row != Col)
                return false;

            for (int i = 0; i < Row; i++)
                for (int j = i + 1; j < Col; j++)
                    if (m_data[i, j] != m_data[j, i])
                        return false;

            return true;
        }

        //一阶矩阵->实数
        public double ToDouble()
        {
            Trace.Assert(Row == 1 && Col == 1);

            return m_data[0, 0];
        }

        //convert to string
        public override string ToString()
        {

            string s = "";
            for (int i = 0; i < Row; i++)
            {
                for (int j = 0; j < Col; j++)
                    s += string.Format("{0} ", m_data[i, j]);

                s += "\r\n";
            }
            return s;

        }
       
        //判断矩阵中所有元素是否大约0
        public bool AllValueAboveZeros()
        {
            bool tempbool = false;
            for (int i=0; i< Row;i++)
            {
                for(int j=0; j<Col;j++)
                    if (m_data[i, j] < 0)
                    {
                        tempbool = false;
                        break;
                    }
                    else
                    {
                        tempbool = true;
                    }
            }
            return tempbool;
        }

        //合并矩阵,按列合并矩阵；
        public Matrix CombineMatrixCol(Matrix A)
        {
            if (A.Row != Row)
            {
                System.Exception e = new Exception("两个矩阵行数不相等；");
                throw e;
            }

            Matrix C = new Matrix(Row, A.Col + Col);

           //double[,] m_data = new double[A.Row, A.Col + Col];
            for (int i = 0; i < Row; i++)
             {
                for (int j = 0; j < Col; j++)
                {
                    C.SetMatrixValue(i, j, m_data[i, j]);
                }
                C.SetMatrixRow(i, Col, A.Col + Col - 1,i,0,A.Col-1,A);
             }
            return C;
        }

        //按行合并矩阵
        public Matrix CombineMatrixRow(Matrix A)
        {
            if (A.Col != Col)
            {
                System.Exception e = new Exception("两个矩阵列数不相等；");
                throw e;
            }

            Matrix C = new Matrix(A.Row+Row, Col);

            //double[,] m_data = new double[A.Row, A.Col + Col];
            for (int i = 0; i < Row; i++)
            {
                for (int j = 0; j < Col; j++)
                {
                    C.SetMatrixValue(i, j, m_data[i, j]);
                }
            }
            for (int i = 0; i < A.Row; i++)
            {
                C.SetMatrixRow(i + Row, 0, Col - 1, i, 0, Col - 1, A);
            }
               
            return C;
        }


        //Get the minimum value in the specified col;
        public double MinMatrixCol(int col)
        {
            double temp = m_data[0,col];
            if(double.IsNaN(temp))
            {
                temp = double.PositiveInfinity;
            } 
            for(int i = 1;i< Row;i++)
                if (m_data[i, col] < temp)
                {
                    temp = m_data[i, col];
                }
            return temp;
        }
      
        //Get the maximum value in the specified col;
        public double MaxMatrixCol(int col)
        {
            double temp = m_data[0, col];
            if (double.IsNaN(temp))
            {
                temp = double.NegativeInfinity;
            } 
            for (int i = 1; i < Row; i++)
                if (m_data[i, col] > temp)
                {
                    temp = m_data[i, col];
                }
            return temp;
        }

        //Get the minimum value in the specified Row;
        public double MinMatrixRow(int row)
        {
            double temp = m_data[row, 0];
            if (double.IsNaN(temp))
            {
                temp = double.PositiveInfinity;
            } 
            for (int i = 1; i < Col; i++)
                if (m_data[row,i] < temp)
                {
                    temp = m_data[row,i];
                }
            return temp;
        }
        //Get the maximum value in the specified Row;
        public double MaxMatrixRow(int row)
        {
            double temp = m_data[row, 0];
            if (double.IsNaN(temp))
            {
                temp = double.NegativeInfinity;
            } 
            for (int i = 1; i < Col; i++)
                if (m_data[row, i] > temp)
                {
                    temp = m_data[row, i];
                }
            return temp;
        }

        //Get the index of the smaller value between the two col values;
        public int MinIndex(int RowIndex1, int RowIndex2, int ColMatrix)
        {
            if (m_data[RowIndex1, ColMatrix] > m_data[RowIndex2, ColMatrix])
                return RowIndex2;
            else
                return RowIndex1;
        }

        //Get the index of the smallest value in the specified col ;
        public int ColMinIndex(int ColMatrix)
        {
            double temp = m_data[0, ColMatrix];
            int tempindex = 0;
            for (int i = 1; i < Row; i++)
            { 
                if(temp > m_data[i,ColMatrix])
                {
                    tempindex = i;
                    temp = m_data[i, ColMatrix];
                }
            }
            return tempindex;
        }


        //Get the row(both singlerow  and matrix) 
        public Matrix GetMatrixRow(int row,int startcol,int endcol)
        {
            Matrix a = new Matrix(endcol - startcol + 1);
            for (int j = 0; j < endcol -startcol + 1; j++)
            {
                a.SetMatrixValue(j, m_data[row, j+startcol]);
            }
            return a;
        }

        //Set the specified row from the one dimensional matrix;
        public void SetMatrixOneRow(int row, int startcol, int endcol, int Astartcol, int Aendcol, Matrix A)
        {
            if (endcol - startcol != Aendcol - Astartcol)
            {
                System.Exception e = new Exception("赋值的行的列数不符合；");
                throw e;
            }

            for (int j = 0; j < endcol - startcol + 1; j++)
            {
                m_data[row, j + startcol] = A.GetMatrixValue(0,Astartcol + j);
            }
        }


        //Set the specified row with the specified row from the multi-matrix
        public void SetMatrixRow(int row,int startcol,int endcol,int Arow,int Astartcol,int Aendcol,Matrix A)
        {
            if (endcol-startcol != Aendcol-Astartcol)
            {
                System.Exception e = new Exception("赋值的行的列数不符合；");
                throw e;
            }

            for(int j = 0;j < endcol - startcol + 1;j++)
            {
                m_data[row,j+startcol] = A.GetMatrixValue(Arow,Astartcol+j);
            }
        }

        //Get the sum of the specified matrix col;
        public double GetColSum(int col)
        {
            double sum = 0;
            for (int i = 0; i < Row; i++)
            {
                sum = sum + m_data[i, col];
            }
            return sum;
        }

        //Get the sum of the specified matrix Row;
        public double GetRowSum(int Row)
        {
            double sum = 0;
            for (int i = 0; i < Col; i++)
            {
                sum = sum + m_data[Row, i];
            }
            return sum;
        }

        //Get the rank of the matrix
        public int rank()
        {
            Matrix A = this.Gauss();
            int rank = 0;

            for (int i = 0; i < Row; i++)
            {
                for (int j = 0; j < Col; j++)
                {
                    if (A.m_data[i, j] != 0)
                    {
                        rank++;
                        break;
                    }
                }
            }
            return rank;
        }

        //Get the square sum of the matrix
        public double SquareSum()
        {
            double squaresum = 0; ;
            for (int i = 0; i < Col; i++)
            {
                for (int j = 0; j < Row; j++)
                {
                    squaresum += Math.Pow(this.m_data[j, i], 2); 
                }
            }
            return squaresum;
        }
        


        //Private data
        private double[,] m_data;

    }
}//end class Matrix



