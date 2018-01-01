using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MatrixHandle;
using System.Windows;
using Copper_OBBS_Simulation_Platform;
using System.Collections; //It contains the definition of the ArrayList;

namespace PSO
{
    public partial class COPSO
    {
        //Construction
        public COPSO()
        {
        }
        public COPSO(InputInfo DivInputPara)
        {
            InputPara = DivInputPara;
        }
       
        //Get the max variable scope by the matrix A and b;
        public Matrix GetVarScope(Matrix A, Matrix b)
        {
            Matrix VarScope = new Matrix(A.Col);
            Matrix B = new Matrix(A.Row, A.Col);
            int Row = A.Row; int Col = A.Col;
            for (int i = 0; i < Row; i++)
            {
                for (int j = 0; j < Col; j++)
                {

                    if ((-esp1 > A.GetMatrixValue(i, j)) || (esp1 < A.GetMatrixValue(i, j)))
                    {
                        B.SetMatrixValue(i, j, b.GetMatrixValue(i, 0) / A.GetMatrixValue(i, j));
                    }
                }
            }

            if (1 == Row)
            {
                for (int j = 0; j < Col; j++)
                {
                    VarScope.SetMatrixValue(j, B.GetMatrixValue(Row-1, j));
                }
            }
            else
            {

                for (int j = 0; j < Col; j++)
                {
                    for (int i = 0; i < Row; i++)
                    {
                        if (-esp1 < B.GetMatrixValue(i, j) && B.GetMatrixValue(i, j) < esp1)
                        {
                            B.SetMatrixValue(i, j, double.NaN);
                        }

                    }
                }

                for (int j = 0; j < Col; j++)
                {
                    VarScope.SetMatrixValue(j, B.MinMatrixCol(j));
                }
            } 
            return VarScope;
        }
       
        //Initiation of the Swarm;
        public Matrix InitSwarm(int ParticleSize, Matrix A,Matrix b)
        {  
         //  Matrix OneDimParSwarm = new Matrix(ParticleSize);
           int a = A.Row;        //当A是形如[1,1,0,0;
                                //          0,0,1,0;
                               //           0,0,0,1]即第（2，2）位置并不是1的时候，情况未考虑；
           Matrix A1 = new Matrix(a,a);  //Matrix A is seperated into A1 and A2;
           Matrix A2 = new Matrix(a,ParticleSize-a);
           for (int i = 0; i < a ; i++)
           {
               for (int j = 0; j < a; j++)
                   A1.SetMatrixValue(i, j, A.GetMatrixValue(i, j));
           }
           for (int i = 0; i < a ; i++)
           {
               for (int j = a; j < ParticleSize; j++)
                   A2.SetMatrixValue(i, j - a, A.GetMatrixValue(i, j));
           }
           
           Matrix VarScope = GetVarScope(A,b);
           Matrix x1 = new Matrix(ParticleSize - a);    //Used to storage the non-base vector in the Ax = b;     
           Matrix x2 = new Matrix(a);                   //Used to storage the base vector in the Ax = b;   
           bool  Flags = true;
           while (Flags)
            { 
                for (int j = 0; j < ParticleSize-a; j++)
                {
                   // Random rd = new Random();
                    double rand = rd.NextDouble();
                    x1.SetMatrixValue(j,rand*0.8*VarScope.GetMatrixValue(j+a));
                }
                Matrix C = A1.InverseMatrix();
                x2 = A1.InverseMatrix()*b -  A1.InverseMatrix()*A2*x1.Transpose();
                x2 = x2.Transpose();
                if (x2.AllValueAboveZeros())
                {
                     //X2 is not permitted to have members that below zero;
                     break;
                }                           
                else
                {
                    Flags=true;
                }
            }
           Matrix x = x2.CombineMatrixCol(x1);
            return x;
         }
        
        //Base step updating of particle swarm;
        public void BaseStepPSO(Matrix VarScope,Matrix ParSwarm, Matrix OptSwarm, Matrix A, Myfunc GibbsFunc, int k, int Iteration)
        { 
            int ParRow = ParSwarm.Row;
            int ParCol = ParSwarm.Col;ParCol = (ParCol-1)/2;
            Matrix SubTract1=new Matrix(ParRow,ParCol);  
            Matrix SubTract2=new Matrix(ParRow,ParCol);
            Matrix SubTract3=new Matrix(ParRow,ParCol);
            //Matrix SubTract4=new Matrix(ParRow,ParCol);
            //Matrix SubTract5=new Matrix(ParRow,ParCol);
            for (int i = 0;i < ParRow; i++)
            {
                //Get the substraction of particle's historical best location and the current location, and storage it in the SubTract1 matrix；
                SubTract1.SetMatrixOneRow(i, 0, ParCol - 1, 0,ParCol-1,OptSwarm.GetMatrixRow(i, 0, ParCol - 1) - ParSwarm.GetMatrixRow(i, 0, ParCol - 1));

                //Get the substraction of swarm's historical best location and the current location, and storage it in the SubTract2 matrix；
                SubTract2.SetMatrixOneRow(i, 0, ParCol - 1,0, ParCol - 1, OptSwarm.GetMatrixRow(ParRow, 0, ParCol - 1) - ParSwarm.GetMatrixRow(i, 0, ParCol - 1));
                
                //Get the Localbest particle in the single-link topology with one neighbourhood and save its substraction in the Subtract3; COPSO uses theselfless model；
                if(i==0)
                {
                  int tempRow = ParSwarm.MinIndex(i+1,ParRow-2,2*ParCol);
                  SubTract3.SetMatrixOneRow(i, 0, ParCol - 1,0, ParCol - 1, ParSwarm.GetMatrixRow(tempRow, 0, ParCol - 1) - ParSwarm.GetMatrixRow(i, 0, ParCol - 1));
                  
                }
                else if(i==1)
                {
                    int tempRow = ParSwarm.MinIndex(i+1,ParRow-1,2*ParCol);
                    SubTract3.SetMatrixOneRow(i, 0, ParCol - 1,0, ParCol - 1, ParSwarm.GetMatrixRow(tempRow, 0, ParCol - 1) - ParSwarm.GetMatrixRow(i, 0, ParCol - 1));
                }
                else if(i==ParRow-1)
                {
                    int tempRow = ParSwarm.MinIndex(0,i-2,2*ParCol);
                    SubTract3.SetMatrixOneRow(i, 0, ParCol - 1,0, ParCol - 1, ParSwarm.GetMatrixRow(tempRow, 0, ParCol - 1) - ParSwarm.GetMatrixRow(i, 0, ParCol - 1));
                }
                else
                {
                    int tempRow = ParSwarm.MinIndex(i+1,i-2,2*ParCol);
                    SubTract3.SetMatrixOneRow(i, 0, ParCol - 1,0, ParCol - 1, ParSwarm.GetMatrixRow(tempRow, 0, ParCol - 1) - ParSwarm.GetMatrixRow(i, 0, ParCol - 1));
                } 
            }

            for(int i = 0;i < ParRow;i++)
            {
                double c1 = rd.NextDouble();double c2 = rd.NextDouble();double c=c1+c2;
                double w= Math.Exp(-30*Math.Pow(k/Iteration,10));
                double s = 2 * rd.NextDouble()/Math.Abs(2.0-c-Math.Sqrt(c * (4 - c)));
                    
                Matrix r = new Matrix(ParCol);
                for(int j = 0;j < ParCol;j++)
                {
                    //Velocity update;
                     if(k < Math.Round(0.9*Iteration))
                     {
                         ParSwarm.SetMatrixValue(i,j+ParCol,s*w*ParSwarm.GetMatrixValue(i,j+ParCol)+s*c1*SubTract1.GetMatrixValue(i,j)+s*c2*SubTract3.GetMatrixValue(i,j));
                     }
                     else
                     {
                         ParSwarm.SetMatrixValue(i,j+ParCol,s*w*ParSwarm.GetMatrixValue(i,j+ParCol)+s*c1*SubTract1.GetMatrixValue(i,j)+s*c2*SubTract2.GetMatrixValue(i,j));
                     }
                    //Velocity confine;
                    if(ParSwarm.GetMatrixValue(i,j+ParCol) >esp2 )
                    {
                        r.SetMatrixValue(j,(VarScope.GetMatrixValue(j)-ParSwarm.GetMatrixValue(i,j))/ParSwarm.GetMatrixValue(i,j+ParCol));
                    }   
                    else if(ParSwarm.GetMatrixValue(i,j+ParCol) <-esp2)
                    {
                        r.SetMatrixValue(j,(-ParSwarm.GetMatrixValue(i,j))/ParSwarm.GetMatrixValue(i,j+ParCol));
                    }
                    else
                    {
                        r.SetMatrixValue(j,0);    
                    }
                }
                double rtemp = r.MinMatrixRow(0);
                if (rtemp<1||rd.NextDouble()<0.5)
                {
                    rtemp=rd.NextDouble()*rtemp;
                }
                else
                {
                    rtemp=1;
                }

                Matrix VelocityTemp = rtemp * ParSwarm.GetMatrixRow(i, ParCol, 2 * ParCol - 1);
                
                ParSwarm.SetMatrixOneRow(i, ParCol, 2 * ParCol - 1,0, ParCol - 1, VelocityTemp);  //It satisfies Av = 0 after the velocity updating;
                ParSwarm.SetMatrixOneRow(i, 0, ParCol - 1,0,ParCol-1, ParSwarm.GetMatrixRow(i, 0, ParCol - 1) + ParSwarm.GetMatrixRow(i, ParCol, 2 * ParCol - 1));  //Location update；
                ParSwarm.SetMatrixValue(i, 2 * ParCol, GibbsFunc(ParSwarm.GetMatrixRow(i, 0, ParCol - 1),InputPara));//Get the adaptive value；            
                if (ParSwarm.GetMatrixValue(i, 2 * ParCol) < GibbsFunc(OptSwarm.GetMatrixRow(i, 0, ParCol - 1),InputPara))
                {
                    OptSwarm.SetMatrixOneRow(i, 0, ParCol - 1, 0, ParCol - 1, ParSwarm.GetMatrixRow(i, 0, ParCol - 1));   //If the adaptive value is better, update the OptSwarm matrix；  
                } 
            }

            //Swarm best particle's location update；
            int MinRow=ParSwarm.ColMinIndex(2*ParCol);
            if (GibbsFunc(ParSwarm.GetMatrixRow(MinRow, 0, ParCol - 1),InputPara) < GibbsFunc(OptSwarm.GetMatrixRow(ParRow, 0, ParCol - 1),InputPara))
             {
                 OptSwarm.SetMatrixOneRow(ParRow, 0, ParCol - 1,0, ParCol - 1, ParSwarm.GetMatrixRow(MinRow, 0, ParCol - 1));   
             }


            // Adding perturbation，C-Perturbation and M-Perturbation; It operates with the swarm pbest values;
            // Setting the probability of the perturbation;
            double p =0;
             //if(k < Math.Round(Iteration*0.85))
             //{
                 p = 1-k/Iteration;
             //}

             #region // DE perturbation, using differential evolving to  make a crossover among PBests;

             Matrix temp1 = new Matrix(ParCol);
             if(rd.NextDouble() < p) 
             {
                for (int i = 0;i < ParRow;i++)
                {
                    double r = rd.NextDouble();
                    int p1=i;
                    Random ri=new Random(i);
                    int p2 = ri.Next(0,i);
                    int p3 = ri.Next(0,i);
                    temp1.SetMatrixOneRow(0, 0, ParCol - 1,0,ParCol-1, OptSwarm.GetMatrixRow(p1, 0, ParCol - 1) +
                        r * (OptSwarm.GetMatrixRow(p2, 0, ParCol - 1) - OptSwarm.GetMatrixRow(p3, 0, ParCol - 1)));
                    if (GibbsFunc(temp1.GetMatrixRow(0, 0, ParCol - 1),InputPara) < GibbsFunc(OptSwarm.GetMatrixRow(i, 0, ParCol - 1),InputPara))
                    {
                        OptSwarm.SetMatrixOneRow(i, 0, ParCol - 1, 0,ParCol - 1, temp1.GetMatrixRow(0, 0, ParCol - 1));
                    }
                }
             }
             #endregion

             #region //v perturbation
             //v perturbation, using the idea that Av=0, to add small v perturbation to the current partilces for the purpose improving the local search capability.
             // this perturbation is begin during the last 10 percent of iteration.
             Matrix r1 = new Matrix(ParCol); 
             if (k < Iteration * 0.5)
             {
                 for(int i = 0; i < ParRow; i++)
                 {  
                     Matrix vtemp = new Matrix(ParCol, ParCol);
                     for (int s = 0; s < ParCol; s++)
                     { 
                        Random ri=new Random(i);
                        int TemoRow = ri.Next(1,ParRow);
                        vtemp.SetMatrixRow(s, 0, ParCol - 1, TemoRow, ParCol, 2 * ParCol - 1, ParSwarm);
                     }
                     
                     Matrix v = ParSwarm.GetMatrixRow(i, ParCol, 2 * ParCol - 1)*vtemp;

                         for (int j = 0; j < ParCol; j++)
                         {
                             if (v.GetMatrixValue(j) > esp2)
                                 r1.SetMatrixValue(j, (VarScope.GetMatrixValue(j) - OptSwarm.GetMatrixValue(i, j)) / v.GetMatrixValue(j));
                             else if (v.GetMatrixValue(j) < -esp2)
                             {
                                 r1.SetMatrixValue(j, (-OptSwarm.GetMatrixValue(i, j)) / v.GetMatrixValue(j));
                             }
                             else
                             {
                                 r1.SetMatrixValue(j, 0);
                             }
                         }
                    double rt = r1.MinMatrixRow(0);
                    if (rt < 1 || rd.NextDouble() < 0.5)
                    {
                        rt = rd.NextDouble() * rt;
                    }
                    Matrix temp2 = OptSwarm.GetMatrixRow(i,0,ParCol-1) + rt*v;
                    if (GibbsFunc(temp2.GetMatrixRow(0, 0, ParCol - 1),InputPara) < GibbsFunc(OptSwarm.GetMatrixRow(i, 0, ParCol - 1),InputPara))
                    {
                        OptSwarm.SetMatrixOneRow(i, 0, ParCol - 1, 0, ParCol - 1, temp2.GetMatrixRow(0, 0, ParCol - 1));   //如果个体较优，则跟新个体历史最优值；  
                    } 
                 }
             }
             #endregion

            // M-pertubation is not achieved yet;   

            //The funuction data in C# is delivered by reference types,so the ParSwarm and OptSwarm is updata on each iteration;
            // And the BaseSetpPSO function do not need return any values;
            
            //Matrix Swarm = new Matrix(ParRow+1,3*ParCol+1);
            //for (int i = 0;i < ParRow;i++)
            //{
            //    Swarm.SetMatrixOneRow(i, 0, 2 * ParCol,0, 2*ParCol, ParSwarm.GetMatrixRow(i, 0, 2 * ParCol));
            //    Swarm.SetMatrixOneRow(i, 2 * ParCol + 1, 3 * ParCol, 0, ParCol - 1, OptSwarm.GetMatrixRow(i, 0, ParCol - 1));
            //}
            //Swarm.SetMatrixOneRow(ParRow, 2 * ParCol + 1, 3 * ParCol, 0, ParCol - 1, OptSwarm.GetMatrixRow(ParRow, 0, ParCol - 1));
           // return Swarm;
        }
      
        //Adapter Value Function;
        public double GibbsFunc(Matrix x, InputInfo InputPara)
        {
            //It should be mentioned that the T1 and T2 is suggested to be double type;
            //If choose to an int type, it maybe ingore the decimal when two int number divided.
            double T1 = 1473;                                               //The Matte and Slag temperature;
            double T2 = 1573;                                               //The Gas temperature;
            double S1 = 8.314 * T1;                                        //S1=R*T1
            double S2 = 8.314 * T2;
            double NCaOs, NMgOs, NAl2O3s, NSiO2s, NO2g, NN2g, NH2Og;
            //InputPara.INPUTMOLE = new double[15]{CuMole, FeMole, SMole, PbMole, ZnMole, AsMole, SbMole, BiMole,
            //                                            CaOMole,MgOMole,Al2O3Mole,SiO2Mole,O2Mole,N2Mole,H2OMole};
            if (InputPara.TEMPINPUTMOLE == null)
            {
                 NCaOs = InputPara.INPUTMOLE[8]; 
                 NMgOs = InputPara.INPUTMOLE[9];
                 NAl2O3s = InputPara.INPUTMOLE[10];
                 NSiO2s = InputPara.INPUTMOLE[11];
                 NO2g = InputPara.O2PERCENT * InputPara.INPUTMOLE[12];
                 NN2g = InputPara.INPUTMOLE[13];
                 NH2Og = InputPara.INPUTMOLE[14];
            }
            else
            {
                NCaOs = InputPara.TEMPINPUTMOLE[8];
                NMgOs = InputPara.TEMPINPUTMOLE[9];
                NAl2O3s = InputPara.TEMPINPUTMOLE[10];
                NSiO2s = InputPara.TEMPINPUTMOLE[11];
                NO2g = InputPara.O2PERCENT * InputPara.TEMPINPUTMOLE[12];
                NN2g = InputPara.TEMPINPUTMOLE[13];
                NH2Og = InputPara.TEMPINPUTMOLE[14];
            }

            //    %%Constraint matrix
            //    %%Cu_COPSO: 1Cu2Sm,2Cum,3Cu2Os,
            //    %%Fe_COPSO: 4FeSm,5FeOs,6Fe3O4s,
            //    %%SO_COPSO: 7SO3g,8S2g,9SO2g,
            //    %%As_COPSO: 10Asm,11As2O3s,12As2g,13AsSg,14AsOg,
            //    %%Sb_COPSO: 15Sbm,16Sb2O3s,17SbOg,18SbSg,
            //    %%Bi_COPSO: 19Bim,20Bi2O3s,21BiOg,
            //    %%Pb_COPSO: 22Pbm,23PbSm,24PbOs,25PbOg, 26PbSg
            //    %%Zn_COPSO: 27ZnSm,28ZnOs,29Zng,30ZnOg.
            //Matte Phase
            double NCu2Sm = 0, NCum = esp1, NFeSm = 0, NAsm = 0, NSbm = 0, 
                   NBim = 0, NPbm = 0, NPbSm = 0, NZnSm = 0;
            //Slag Phase
            double NCu2Os = 0, NFeOs = 0, NFe3O4s = 0, NAs2O3s = 0,
                   NSb2O3s = 0, NBi2O3s = 0, NPbOs = 0, NZnOs = 0;
            //Gas Phase
            double NSO3g = 0, NS2g = 0, NSO2g = 0, NAs2g = 0, 
                   NAsSg = 0, NAsOg = 0, NSbOg = 0, NSbSg = 0,
                   NBi2O3g = 0, NPbOg = 0, NPbSg = 0, NZng = 0, NZnSg = 0;
            
                  //{NCu2Sm, NCum, NCu2Os, NFeSm, NFeOs, NFe3O4s,
                  //NSO3g, NS2g, NSO2g, NAsm, NAs2O3s, NAs2g, NAsSg,
                  //NAsOg, NSbm, NSb2O3s, NSbOg, NSbSg, NBim, NBi2O3s,
                  //NBi2O3g, NPbm, NPbSm, NPbOs, NPbOg, NPbSg, NZnSm, NZnOs, NZng,NZnOg};
            int z = 0;
            if (InputPara.SUBSTANCECHOOSE[0] != 0)  { NCu2Sm = x[0,z++];}
            if (InputPara.SUBSTANCECHOOSE[1] != 0)  { NCum = x[0, z++];}
            if (InputPara.SUBSTANCECHOOSE[2] != 0)  { NCu2Os = x[0, z++];}
            if (InputPara.SUBSTANCECHOOSE[3] != 0)  { NFeSm = x[0, z++];}
            if (InputPara.SUBSTANCECHOOSE[4] != 0)  { NFeOs = x[0, z++];}
            if (InputPara.SUBSTANCECHOOSE[5] != 0)  { NFe3O4s = x[0, z++];}
            if (InputPara.SUBSTANCECHOOSE[6] != 0)  { NSO3g = x[0, z++];}
            if (InputPara.SUBSTANCECHOOSE[7] != 0)  { NS2g = x[0, z++];}
            if (InputPara.SUBSTANCECHOOSE[8] != 0)  { NSO2g = x[0, z++];}
            if (InputPara.SUBSTANCECHOOSE[9] != 0)  { NAsm = x[0, z++];}
            if (InputPara.SUBSTANCECHOOSE[10] != 0) { NAs2O3s = x[0, z++];}
            if (InputPara.SUBSTANCECHOOSE[11] != 0) { NAs2g = x[0, z++];}
            if (InputPara.SUBSTANCECHOOSE[12] != 0) { NAsSg = x[0, z++];}
            if (InputPara.SUBSTANCECHOOSE[13] != 0) { NAsOg = x[0, z++];}
            if (InputPara.SUBSTANCECHOOSE[14] != 0) { NSbm = x[0, z++];}
            if (InputPara.SUBSTANCECHOOSE[15] != 0) { NSb2O3s = x[0, z++];}
            if (InputPara.SUBSTANCECHOOSE[16] != 0) { NSbOg = x[0, z++];}
            if (InputPara.SUBSTANCECHOOSE[17] != 0) { NSbSg = x[0, z++];}
            if (InputPara.SUBSTANCECHOOSE[18] != 0) { NBim = x[0, z++];}
            if (InputPara.SUBSTANCECHOOSE[19] != 0) { NBi2O3s = x[0, z++];}
            if (InputPara.SUBSTANCECHOOSE[20] != 0) { NBi2O3g = x[0, z++];}
            if (InputPara.SUBSTANCECHOOSE[21] != 0) { NPbm = x[0, z++];}
            if (InputPara.SUBSTANCECHOOSE[22] != 0) { NPbSm = x[0, z++];}
            if (InputPara.SUBSTANCECHOOSE[23] != 0) { NPbOs = x[0, z++];}
            if (InputPara.SUBSTANCECHOOSE[24] != 0) { NPbOg = x[0, z++];}
            if (InputPara.SUBSTANCECHOOSE[25] != 0) { NPbSg = x[0, z++]; }
            if (InputPara.SUBSTANCECHOOSE[26] != 0) { NZnSm = x[0, z++];}
            if (InputPara.SUBSTANCECHOOSE[27] != 0) { NZnOs = x[0, z++];}
            if (InputPara.SUBSTANCECHOOSE[28] != 0) { NZng = x[0, z++]; }
            if (InputPara.SUBSTANCECHOOSE[29] !=0)  {NZnSg = x[0,z++];}
            double Cu = (63.5*2*NCu2Sm + 63.5*NCum)/(159*NCu2Sm+ 63.5*NCum+88*NFeSm+74.92*NAsm+121.8*NSbm+
                         209*NBim+207.2*NPbm+237.2*NPbSm+97.39*NZnSm);
                   Cu = Cu*100;
            double NMatte = NCu2Sm + NCum + NFeSm + NAsm + NSbm + NBim + NPbm + NPbSm + NZnSm;
            double NSlag = NCu2Os + NFeOs + NFe3O4s + NAs2O3s + NSb2O3s + NBi2O3s + NPbOs + NZnOs + NCaOs + NMgOs + NAl2O3s + NSiO2s;
            double NGas = NO2g + NSO3g + NS2g + NSO2g + NAs2g + NAsSg + NAsOg + NSbOg + NSbSg + NBi2O3g + NPbOg +NPbSg + NZng + NZnSg + NN2g + NH2Og;
            


            //double yt3 = NCu2Sm * (-122884 + S1 * Math.Log(NCu2Sm / NMatte));
            //double y1 =  NCum * (-71.94 + Math.Log(NCum / NMatte));
            //double yt2 = NFeSm * (-103104 + S1 * Math.Log(NFeSm / NMatte * (0.925 / ((NCu2Sm / NMatte) + 1))));
            //double y2 = NFeOs * (-167720 + S1 * Math.Log(NFeOs / NSlag * (1.42 * (NFeOs / NSlag) - 0.044)));
            //double y3 = 2.25 * NFe3O4s * (-648236 + S1 * Math.Log(NFe3O4s / NSlag * (0.69 + 56.8 * (NFe3O4s / NSlag) + 5.45 * (NSiO2s / NSlag))));
            //double y4 = NSO2g * (-2863291 + S2 * Math.Log(NSO2g / NGas));

            double y = NCu2Sm * (-122885 + S1 * Math.Log(NCu2Sm / NMatte)) +
                       NCum * (-57.10 + Math.Log(NCum / NMatte * 14)) +
                       NFeSm * (-103104 + S1 * Math.Log(NFeSm / NMatte * (0.925 / ((NCu2Sm / NMatte) + 1)))) +
                       //NFeOm * (-167720 + S1 * Math.Log(NFeOm / NMatte * (Math.Exp(5.1 + 6.2 * Math.Log(NCu2Sm / NMatte) + 6.41 * Math.Pow((Math.Log(NCu2Sm / NMatte)) , 2) + 2.8 *Math.Pow((Math.Log(NCu2Sm / NMatte)), 3))))) +
                       // 2.25*NFe3O4m * (-648236 + S1 * Math.Log(NFe3O4m / NMatte * (Math.Exp(4.96 + 9.9 * Math.Log(NCu2Sm / NMatte) + 7.43 *Math.Pow((Math.Log(NCu2Sm / NMatte)) , 2) + 2.55 * Math.Pow((Math.Log(NCu2Sm / NMatte)) , 3))))) +
                       NFeOs * (-167720 + S1 * Math.Log(NFeOs / NSlag * (1.42 * (NFeOs / NSlag) - 0.044))) +
                       NFe3O4s * (-648236 + S1 * Math.Log(NFe3O4s / NSlag * (0.69 + 56.8 * (NFe3O4s / NSlag) + 5.45 * (NSiO2s / NSlag)))) +
                       NSiO2s * (-649130 + S1 * Math.Log(NSiO2s / NSlag * 2.1)) +
                       //NCu2Ss * (-122884 + S1 * Math.Log(NCu2Ss / NSlag * (Math.Exp(2.46 + 6.22 * (NCu2Sm / NMatte))))) +
                       1.13*NCu2Os * (-59872 + S1 * Math.Log(NCu2Os / NSlag * (57.14 * NCu2Os / NSlag))) +
                       //NFeSs*(-103104+S1 *Math.Log(NFeSs/NSlag*70))+
                       NSO2g * (-286329 + S2 * Math.Log(NSO2g / NGas)) +
                       NS2g * (-77528 + S2 * Math.Log(NS2g / NGas)) +
                       NO2g * (S2 * Math.Log(NO2g / NGas)) +
                       NSO3g * (-238954 + S2 * Math.Log(NSO3g / NGas)) +
                       0.63*NAsm * (S1 * Math.Log(NAsm / NMatte * (8.087 - 0.128 * Cu + 0.014 * Cu * Math.Log10(Cu)))) +
                       1.127*NAs2O3s * (-318303 + S1 * Math.Log(NAs2O3s / NSlag * 3.838 * Math.Exp(1523 / T1) * Math.Pow((NO2g / NGas), 0.158))) +
                       NAs2g * (-32184 + S2 * Math.Log(NAs2g / NGas)) +
                       NAsSg * (-48779 + S2 * Math.Log(NAsSg / NGas)) +
                       NAsOg * (-82957 + S2 * Math.Log(NAsOg / NGas)) +
                       3.343*NSbm * (S1 * Math.Log(NSbm / NMatte * (-0.996 + 2.42 * Cu - 1.26 * Cu * Math.Log10(Cu)))) +
                       1.768*NSb2O3s * (-356262 + S1 * Math.Log(NSb2O3s / NSlag * Math.Exp(1055.66 / T1))) +
                       NSbOg * (-220636 + S2 * Math.Log(NSbOg / NGas)) +
                       NSbSg * (-28843 + S2 * Math.Log(NSbSg / NGas)) +
                       1.069*NBim * (S1 * Math.Log(NBim / NMatte * Math.Pow(10, (1900 / T1 - 0.464)))) +
                       // NBim * (S1 * Math.Log(NBim / NMatte * Math.Exp(-3.657+9337/T1-(0.863-1066/T1)*Cu) ))+
                       1.497*NBi2O3s * (-187895 + S1 * Math.Log(NBi2O3s / NSlag * Math.Exp(-1055.66 / T1))) +
                       2.98*NBi2O3g * (-53168 + S2 * Math.Log(NBi2O3g / NGas)) +
                       1.1148*NPbm * (S1 * Math.Log(NPbm / NMatte * 23)) +
                       1.1148*NPbSm * (-70326 + S1 * Math.Log(NPbSm / NMatte * Math.Exp((-2.716 + 2441 / T1 + (0.815 - 3610 / T1) * (80 - Cu) / 100)))) +
                       0.995*NPbOs * (-80851 + S1 * Math.Log(NPbOs / NSlag * Math.Exp(-3330 / T1))) +
                       NPbOg * (-23612 + S2 * Math.Log(NPbOg / NGas)) +
                       NPbSg *(-64462.759 +S2*Math.Log(NPbSg /NGas))+
                       1.05*NZnSm * (-154702 + S1 * Math.Log(NZnSm / NMatte * Math.Exp((-2.054 + 6917 / T1 - (1.522 - 1032 / T1) * (80 - Cu) / 100)))) +
                       0.89*NZnOs * (-196868 + S1 * Math.Log(NZnOs / NSlag * Math.Exp(920 / T1))) +
                       NZng * (-37941 + S2 * Math.Log(NZng / NGas))+
                       NZnSg*(-13952 +S2*Math.Log(NZnSg/NGas));           
                       // NZnOg * (1144 + S2 * Math.Log(NZnOg / NGas));
                return y;
        }
        
        //public double GibbsFunc(Matrix x, InputInfo InputPara, double[] Percent)
        //{
        //    //It should be mentioned that the T1 and T2 is suggested to be double type;
        //    //If choose to an int type, it maybe ingore the decimal when two int number divided.
        //    double T1 = 1473;                                               //The Matte and Slag temperature;
        //    double T2 = 1573;                                               //The Gas temperature;
        //    double S1 = 8.314 * T1;                                        //S1=R*T1
        //    double S2 = 8.314 * T2;

        //    double magnitude = Math.Pow(10, 4);
        //    //InputPara.INPUTMOLE = new double[15]{CuMole, FeMole, SMole, PbMole, ZnMole, AsMole, SbMole, BiMole,
        //    //                                            CaOMole,MgOMole,Al2O3Mole,SiO2Mole,O2Mole,N2Mole,H2OMole};
        //    double NCaOs = Percent[8] * Percent[13] / 56 * magnitude + Percent[15] * 2 / (40 + 16) * magnitude;
        //    double NMgOs = Percent[9] * Percent[13] / 40.3 * magnitude + Percent[15] * 2 / (24.3 + 12) * magnitude;
        //    double NAl2O3s = Percent[10] * Percent[13] / (27 * 2 + 16 * 3) * magnitude + Percent[15] * 1 / (27 * 2 + 16 * 3) * magnitude;
        //    double  NSiO2s = Percent[11] * Percent[13] / 60 * magnitude + Percent[15] * 95 / (28.08 + 32) * magnitude;
        //    double NO2g = InputPara.O2PERCENT * (Percent[16] + Percent[17] * 0.21) * 1.4276 / 32 * 1000 * 2;
        //    double NN2g = Percent[17] * 0.79 * 1.24988 / 28 * 1000;
        //    double NH2Og = Percent[13] * Percent[14] / (100 - Percent[14]) / 18 * 100 * magnitude; 
        //    //    %%Constraint matrix
        //    //    %%Cu_COPSO: 1Cu2Sm,2Cum,3Cu2Os,
        //    //    %%Fe_COPSO: 4FeSm,5FeOs,6Fe3O4s,
        //    //    %%SO_COPSO: 7SO3g,8S2g,9SO2g,
        //    //    %%As_COPSO: 10Asm,11As2O3s,12As2g,13AsSg,14AsOg,
        //    //    %%Sb_COPSO: 15Sbm,16Sb2O3s,17SbOg,18SbSg,
        //    //    %%Bi_COPSO: 19Bim,20Bi2O3s,21BiOg,
        //    //    %%Pb_COPSO: 22Pbm,23PbSm,24PbOs,25PbOg, 26PbSg
        //    //    %%Zn_COPSO: 27ZnSm,28ZnOs,29Zng,30ZnOg.
        //    //Matte Phase
        //    double NCu2Sm = 0, NCum = esp1, NFeSm = 0, NAsm = 0, NSbm = 0,
        //           NBim = 0, NPbm = 0, NPbSm = 0, NZnSm = 0;
        //    //Slag Phase
        //    double NCu2Os = 0, NFeOs = 0, NFe3O4s = 0, NAs2O3s = 0,
        //           NSb2O3s = 0, NBi2O3s = 0, NPbOs = 0, NZnOs = 0;
        //    //Gas Phase
        //    double NSO3g = 0, NS2g = 0, NSO2g = 0, NAs2g = 0,
        //           NAsSg = 0, NAsOg = 0, NSbOg = 0, NSbSg = 0,
        //           NBiOg = 0, NPbOg = 0, NPbSg = 0, NZng = 0, NZnSg = 0;

        //    //{NCu2Sm, NCum, NCu2Os, NFeSm, NFeOs, NFe3O4s,
        //    //NSO3g, NS2g, NSO2g, NAsm, NAs2O3s, NAs2g, NAsSg,
        //    //NAsOg, NSbm, NSb2O3s, NSbOg, NSbSg, NBim, NBi2O3s,
        //    //NBiOg, NPbm, NPbSm, NPbOs, NPbOg, NPbSg, NZnSm, NZnOs, NZng,NZnOg};
        //    int z = 0;
        //    if (InputPara.SUBSTANCECHOOSE[0] != 0) { NCu2Sm = x[0, z++]; }
        //    if (InputPara.SUBSTANCECHOOSE[1] != 0) { NCum = x[0, z++]; }
        //    if (InputPara.SUBSTANCECHOOSE[2] != 0) { NCu2Os = x[0, z++]; }
        //    if (InputPara.SUBSTANCECHOOSE[3] != 0) { NFeSm = x[0, z++]; }
        //    if (InputPara.SUBSTANCECHOOSE[4] != 0) { NFeOs = x[0, z++]; }
        //    if (InputPara.SUBSTANCECHOOSE[5] != 0) { NFe3O4s = x[0, z++]; }
        //    if (InputPara.SUBSTANCECHOOSE[6] != 0) { NSO3g = x[0, z++]; }
        //    if (InputPara.SUBSTANCECHOOSE[7] != 0) { NS2g = x[0, z++]; }
        //    if (InputPara.SUBSTANCECHOOSE[8] != 0) { NSO2g = x[0, z++]; }
        //    if (InputPara.SUBSTANCECHOOSE[9] != 0) { NAsm = x[0, z++]; }
        //    if (InputPara.SUBSTANCECHOOSE[10] != 0) { NAs2O3s = x[0, z++]; }
        //    if (InputPara.SUBSTANCECHOOSE[11] != 0) { NAs2g = x[0, z++]; }
        //    if (InputPara.SUBSTANCECHOOSE[12] != 0) { NAsSg = x[0, z++]; }
        //    if (InputPara.SUBSTANCECHOOSE[13] != 0) { NAsOg = x[0, z++]; }
        //    if (InputPara.SUBSTANCECHOOSE[14] != 0) { NSbm = x[0, z++]; }
        //    if (InputPara.SUBSTANCECHOOSE[15] != 0) { NSb2O3s = x[0, z++]; }
        //    if (InputPara.SUBSTANCECHOOSE[16] != 0) { NSbOg = x[0, z++]; }
        //    if (InputPara.SUBSTANCECHOOSE[17] != 0) { NSbSg = x[0, z++]; }
        //    if (InputPara.SUBSTANCECHOOSE[18] != 0) { NBim = x[0, z++]; }
        //    if (InputPara.SUBSTANCECHOOSE[19] != 0) { NBi2O3s = x[0, z++]; }
        //    if (InputPara.SUBSTANCECHOOSE[20] != 0) { NBiOg = x[0, z++]; }
        //    if (InputPara.SUBSTANCECHOOSE[21] != 0) { NPbm = x[0, z++]; }
        //    if (InputPara.SUBSTANCECHOOSE[22] != 0) { NPbSm = x[0, z++]; }
        //    if (InputPara.SUBSTANCECHOOSE[23] != 0) { NPbOs = x[0, z++]; }
        //    if (InputPara.SUBSTANCECHOOSE[24] != 0) { NPbOg = x[0, z++]; }
        //    if (InputPara.SUBSTANCECHOOSE[25] != 0) { NPbSg = x[0, z++]; }
        //    if (InputPara.SUBSTANCECHOOSE[26] != 0) { NZnSm = x[0, z++]; }
        //    if (InputPara.SUBSTANCECHOOSE[27] != 0) { NZnOs = x[0, z++]; }
        //    if (InputPara.SUBSTANCECHOOSE[28] != 0) { NZng = x[0, z++]; }
        //    if (InputPara.SUBSTANCECHOOSE[29] != 0) { NZnSg = x[0, z++]; }
        //    double Cu = (63.5 * 2 * NCu2Sm + 63.5 * NCum) / (159 * NCu2Sm + 63.5 * NCum + 88 * NFeSm + 74.92 * NAsm + 121.8 * NSbm +
        //                 209 * NBim + 207.2 * NPbm + 237.2 * NPbSm + 97.39 * NZnSm);
        //    Cu = Cu * 100;
        //    double NMatte = NCu2Sm + NCum + NFeSm + NAsm + NSbm + NBim + NPbm + NPbSm + NZnSm;
        //    double NSlag = NCu2Os + NFeOs + NFe3O4s + NAs2O3s + NSb2O3s + NBi2O3s + NPbOs + NZnOs + NCaOs + NMgOs + NAl2O3s + NSiO2s;
        //    double NGas = NO2g + NSO3g + NS2g + NSO2g + NAs2g + NAsSg + NAsOg + NSbOg + NSbSg + NBiOg + NPbOg + NPbSg + NZng + NZnSg + NN2g + NH2Og;



        //    //double yt3 = NCu2Sm * (-122884 + S1 * Math.Log(NCu2Sm / NMatte));
        //    //double y1 =  NCum * (-71.94 + Math.Log(NCum / NMatte));
        //    //double yt2 = NFeSm * (-103104 + S1 * Math.Log(NFeSm / NMatte * (0.925 / ((NCu2Sm / NMatte) + 1))));
        //    //double y2 = NFeOs * (-167720 + S1 * Math.Log(NFeOs / NSlag * (1.42 * (NFeOs / NSlag) - 0.044)));
        //    //double y3 = 2.25 * NFe3O4s * (-648236 + S1 * Math.Log(NFe3O4s / NSlag * (0.69 + 56.8 * (NFe3O4s / NSlag) + 5.45 * (NSiO2s / NSlag))));
        //    //double y4 = NSO2g * (-2863291 + S2 * Math.Log(NSO2g / NGas));

        //    double y = NCu2Sm * (-122885 + S1 * Math.Log(NCu2Sm / NMatte)) +
        //               NCum * (-57.10 + Math.Log(NCum / NMatte * 14)) +
        //               NFeSm * (-103104 + S1 * Math.Log(NFeSm / NMatte * (0.925 / ((NCu2Sm / NMatte) + 1)))) +
        //        //NFeOm * (-167720 + S1 * Math.Log(NFeOm / NMatte * (Math.Exp(5.1 + 6.2 * Math.Log(NCu2Sm / NMatte) + 6.41 * Math.Pow((Math.Log(NCu2Sm / NMatte)) , 2) + 2.8 *Math.Pow((Math.Log(NCu2Sm / NMatte)), 3))))) +
        //        // 2.25*NFe3O4m * (-648236 + S1 * Math.Log(NFe3O4m / NMatte * (Math.Exp(4.96 + 9.9 * Math.Log(NCu2Sm / NMatte) + 7.43 *Math.Pow((Math.Log(NCu2Sm / NMatte)) , 2) + 2.55 * Math.Pow((Math.Log(NCu2Sm / NMatte)) , 3))))) +
        //               NFeOs * (-167720 + S1 * Math.Log(NFeOs / NSlag * (1.42 * (NFeOs / NSlag) - 0.044))) +
        //               NFe3O4s * (-648236 + S1 * Math.Log(NFe3O4s / NSlag * (0.69 + 56.8 * (NFe3O4s / NSlag) + 5.45 * (NSiO2s / NSlag)))) +
        //               NSiO2s * (-649130 + S1 * Math.Log(NSiO2s / NSlag * 2.1)) +
        //        //NCu2Ss * (-122884 + S1 * Math.Log(NCu2Ss / NSlag * (Math.Exp(2.46 + 6.22 * (NCu2Sm / NMatte))))) +
        //               1.05 * NCu2Os * (-59872 + S1 * Math.Log(NCu2Os / NSlag * (57.14 * NCu2Os / NSlag))) +
        //        //NFeSs*(-103104+S1 *Math.Log(NFeSs/NSlag*70))+
        //               NSO2g * (-286329 + S2 * Math.Log(NSO2g / NGas)) +
        //               NS2g * (-77528 + S2 * Math.Log(NS2g / NGas)) +
        //               NO2g * (S2 * Math.Log(NO2g / NGas)) +
        //               NSO3g * (-238954 + S2 * Math.Log(NSO3g / NGas)) +
        //               1.28 * NAsm * (S1 * Math.Log(NAsm / NMatte * (8.087 - 0.128 * Cu + 0.014 * Cu * Math.Log10(Cu)))) +
        //               1.249 * NAs2O3s * (-318303 + S1 * Math.Log(NAs2O3s / NSlag * 3.838 * Math.Exp(1523 / T1) * Math.Pow((NO2g / NGas), 0.158))) +
        //               NAs2g * (-32184 + S2 * Math.Log(NAs2g / NGas)) +
        //               NAsSg * (-48779 + S2 * Math.Log(NAsSg / NGas)) +
        //               NAsOg * (-82957 + S2 * Math.Log(NAsOg / NGas)) +
        //               4.68 * NSbm * (S1 * Math.Log(NSbm / NMatte * (-0.996 + 2.42 * Cu - 1.26 * Cu * Math.Log10(Cu)))) +
        //               1.706 * NSb2O3s * (-356262 + S1 * Math.Log(NSb2O3s / NSlag * Math.Exp(1055.66 / T1))) +
        //               NSbOg * (-220636 + S2 * Math.Log(NSbOg / NGas)) +
        //               NSbSg * (-28843 + S2 * Math.Log(NSbSg / NGas)) +
        //               2.64 * NBim * (S1 * Math.Log(NBim / NMatte * Math.Pow(10, (1900 / T1 - 0.464)))) +
        //               2.15 * NBi2O3s * (-187895 + S1 * Math.Log(NBi2O3s / NSlag * Math.Exp(-1055.66 / T1))) +
        //               3.23 * NBiOg * (21125 + S2 * Math.Log(NBiOg / NGas)) +
        //               NPbm * (S1 * Math.Log(NPbm / NMatte * 23)) +
        //               NPbSm * (-70326 + S1 * Math.Log(NPbSm / NMatte * Math.Exp((-2.716 + 2441 / T1 + (0.815 - 3610 / T1) * (80 - Cu) / 100)))) +
        //               1.105 * NPbOs * (-80851 + S1 * Math.Log(NPbOs / NSlag * Math.Exp(-3330 / T1))) +
        //               NPbOg * (-23612 + S2 * Math.Log(NPbOg / NGas)) +
        //               1.72 * NPbSg * (S2 * Math.Log(NPbSg / NGas)) +
        //               3.28 * NZnSm * (-154702 + S1 * Math.Log(NZnSm / NMatte * Math.Exp((-2.054 + 6917 / T1 - (1.522 - 1032 / T1) * (80 - Cu) / 100)))) +
        //               2.58 * NZnOs * (-196868 + S1 * Math.Log(NZnOs / NSlag * Math.Exp(920 / T1))) +
        //               3.85 * NZng * (-37941 + S2 * Math.Log(NZng / NGas)) +
        //               3.85 * NZnSg * (-13952 + S2 * Math.Log(NZnSg / NGas));
        //    // NZnOg * (1144 + S2 * Math.Log(NZnOg / NGas));
        //    return y;
        //}
       
        //Parameter of postback;

        public float[,] GetProcessResults(bool IsVarResults, double[] x, InputInfo InputPara, bool IsFeSiO2)
        {
            double NCaOs, NMgOs, NAl2O3s, NSiO2s, NO2g, NN2g, NH2Og;
            if (!IsVarResults)
            {
                NCaOs = InputPara.INPUTMOLE[8];
                NMgOs = InputPara.INPUTMOLE[9];
                NAl2O3s = InputPara.INPUTMOLE[10];
                NSiO2s = InputPara.INPUTMOLE[11];
                NO2g = InputPara.O2PERCENT * InputPara.INPUTMOLE[12];
                NN2g = InputPara.INPUTMOLE[13];
                NH2Og = InputPara.INPUTMOLE[14];
            }
            else 
            {
                NCaOs = InputPara.TEMPINPUTMOLE[8];
                NMgOs = InputPara.TEMPINPUTMOLE[9];
                NAl2O3s = InputPara.TEMPINPUTMOLE[10];
                NSiO2s = InputPara.TEMPINPUTMOLE[11];
                NO2g = InputPara.O2PERCENT * InputPara.TEMPINPUTMOLE[12];
                NN2g = InputPara.TEMPINPUTMOLE[13];
                NH2Og = InputPara.TEMPINPUTMOLE[14];
            }
            
            //Matte Phase
            double NCu2Sm = 0, NCum = 0, NFeSm = 0, NAsm = 0, NSbm = 0,
                   NBim = 0, NPbm = 0, NPbSm = 0, NZnSm = 0;
            //Slag Phase
            double NCu2Os = 0, NFeOs = 0, NFe3O4s = 0, NAs2O3s = 0,
                   NSb2O3s = 0, NBi2O3s = 0, NPbOs = 0, NZnOs = 0;
            //Gas Phase
            double NSO3g = 0, NS2g = 0, NSO2g = 0, NAs2g = 0,
                   NAsSg = 0, NAsOg = 0, NSbOg = 0, NSbSg = 0,
                   NBi2O3g = 0, NPbOg = 0, NPbSg = 0, NZng = 0, NZnOg = 0;
            int z = 0;
            if (InputPara.SUBSTANCECHOOSE[0] != 0) { NCu2Sm = x[z++]; }
            if (InputPara.SUBSTANCECHOOSE[1] != 0) { NCum = x[z++]; }
            if (InputPara.SUBSTANCECHOOSE[2] != 0) { NCu2Os = x[z++]; }
            if (InputPara.SUBSTANCECHOOSE[3] != 0) { NFeSm = x[z++]; }
            if (InputPara.SUBSTANCECHOOSE[4] != 0) { NFeOs = x[z++]; }
            if (InputPara.SUBSTANCECHOOSE[5] != 0) { NFe3O4s = x[z++]; }
            if (InputPara.SUBSTANCECHOOSE[6] != 0) { NSO3g = x[z++]; }
            if (InputPara.SUBSTANCECHOOSE[7] != 0) { NS2g = x[z++]; }
            if (InputPara.SUBSTANCECHOOSE[8] != 0) { NSO2g = x[z++]; }
            if (InputPara.SUBSTANCECHOOSE[9] != 0) { NAsm = x[z++]; }
            if (InputPara.SUBSTANCECHOOSE[10] != 0) { NAs2O3s = x[z++]; }
            if (InputPara.SUBSTANCECHOOSE[11] != 0) { NAs2g = x[z++]; }
            if (InputPara.SUBSTANCECHOOSE[12] != 0) { NAsSg = x[z++]; }
            if (InputPara.SUBSTANCECHOOSE[13] != 0) { NAsOg = x[z++]; }
            if (InputPara.SUBSTANCECHOOSE[14] != 0) { NSbm = x[z++]; }
            if (InputPara.SUBSTANCECHOOSE[15] != 0) { NSb2O3s = x[z++]; }
            if (InputPara.SUBSTANCECHOOSE[16] != 0) { NSbOg = x[z++]; }
            if (InputPara.SUBSTANCECHOOSE[17] != 0) { NSbSg = x[z++]; }
            if (InputPara.SUBSTANCECHOOSE[18] != 0) { NBim = x[z++]; }
            if (InputPara.SUBSTANCECHOOSE[19] != 0) { NBi2O3s = x[z++]; }
            if (InputPara.SUBSTANCECHOOSE[20] != 0) { NBi2O3g = x[z++]; }
            if (InputPara.SUBSTANCECHOOSE[21] != 0) { NPbm = x[z++]; }
            if (InputPara.SUBSTANCECHOOSE[22] != 0) { NPbSm = x[z++]; }
            if (InputPara.SUBSTANCECHOOSE[23] != 0) { NPbOs = x[z++]; }
            if (InputPara.SUBSTANCECHOOSE[24] != 0) { NPbOg = x[z++]; }
            if (InputPara.SUBSTANCECHOOSE[25] != 0) { NPbSg = x[z++]; }
            if (InputPara.SUBSTANCECHOOSE[26] != 0) { NZnSm = x[z++]; }
            if (InputPara.SUBSTANCECHOOSE[27] != 0) { NZnOs = x[z++]; }
            if (InputPara.SUBSTANCECHOOSE[28] != 0) { NZng = x[z++]; }
            if (InputPara.SUBSTANCECHOOSE[29] != 0) { NZnOg = x[z++]; }
            
            //Mid Variable for the solution of the Mechanical entrainment coefficient.
            double MMatte = (159 * NCu2Sm + 63.5 * NCum + 88 * NFeSm + 74.92 * NAsm + 121.8 * NSbm + 209 * NBim + 207.2 * NPbm + 237.2 * NPbSm + 97.39 * NZnSm);
            double TempMatte = (63.5 * 2 * NCu2Sm + 63.5 * NCum) / (159 * NCu2Sm + 63.5 * NCum) * 100;
            double Temp = (63.5 * 2 * NCu2Sm + 63.5 * NCum);
            double Temp1 = 88 * NFeSm + 74.92 * NAsm + 121.8 * NSbm + 209 * NBim + 207.2 * NPbm + 237.2 * NPbSm + 97.39 * NZnSm;
            double Temp1Pct = Temp1 / MMatte*100;
            double t = (159 * NCu2Sm + 63.5 * NCum) / MMatte * 100;
            double Cu = (63.5 * 2 * NCu2Sm + 63.5 * NCum) / MMatte* 100;
            double MSlag = (63.5 * 2 + 16) * NCu2Os + (56 + 16) * NFeOs + (3 * 56 + 16 * 4) * NFe3O4s +
                          (74.9 * 2 + 16 * 3) * NAs2O3s + (121.7 * 2 + 16 * 3) * NSb2O3s + (208.9 * 2 + 16 * 3) * NBi2O3s + (207 + 16) * NPbOs + (65.5 + 16) * NZnOs +
                           NCaOs * (40 + 16) + NMgOs * (24 + 16) + NAl2O3s * (27 * 2 + 16 * 3) + NSiO2s * (28 + 2 * 16);
            double Fe3O4Slag = (56 * 3 + 16 * 4) * NFe3O4s / MSlag * 100;
            #region //Get the Mechanical entrainment coefficient;
            double MatteToSlagCof, SlagToMatteCof;
            MatteToSlagCof = (2.38253 + 0.35747 * Math.Exp((Cu - 50) / 18.20478)) / 100;
            SlagToMatteCof = (1.51891 + (3.10731 - 1.31891) / (1.0 + Math.Pow(10, (54.22768 - Cu) * 0.05103))) / 100;
            #endregion
            if (IsFeSiO2)
            {
                MatteToSlagCof = MatteToSlagCof + (Fe3O4Slag - 25) / 1000;
                SlagToMatteCof = SlagToMatteCof + (Fe3O4Slag - 25) / 1000;
            }
           


             MSlag = (63.5 * 2 + 16) * NCu2Os + (56 + 16) * NFeOs + (3 * 56 + 16 * 4) * NFe3O4s +
                           (74.9 * 2 + 16 * 3) * NAs2O3s + (121.7 * 2 + 16 * 3) * NSb2O3s + (208.9 * 2 + 16 * 3) * NBi2O3s + (207 + 16) * NPbOs + (65.5 + 16) * NZnOs +
                            NCaOs * (40 + 16) + NMgOs * (24 + 16) + NAl2O3s * (27 * 2 + 16 * 3) + NSiO2s * (28 + 2 * 16);

            //MatteToSlagCof : S1; SlagToMatteCof : S2;
            double MMatteInSlag = ( MatteToSlagCof * SlagToMatteCof  * MMatte + MatteToSlagCof * SlagToMatteCof * MSlag - MSlag * MatteToSlagCof ) / (MatteToSlagCof + SlagToMatteCof - 1);
            double MSlagInMatte = ( MatteToSlagCof * SlagToMatteCof  * MMatte + MatteToSlagCof * SlagToMatteCof * MSlag - MMatte * SlagToMatteCof) / (MatteToSlagCof + SlagToMatteCof - 1);
            double MAppSlag = MSlag - MSlagInMatte;
            double MAppMatte = MMatte - MMatteInSlag;

            double NCu2Ss = MMatteInSlag * NCu2Sm / MMatte, NCus = MMatteInSlag * NCum / MMatte, NFeSs = MMatteInSlag * NFeSm / MMatte;
            double NFeOm = MSlagInMatte * NFeOs / MSlag, NFe3O4m = MSlagInMatte * NFe3O4s / MSlag, NSiO2m = MSlagInMatte * NSiO2s / MSlag;
            NCu2Sm -= NCu2Ss; NCum -= NCus; NFeSm -= NFeSs; NFeOs -= NFeOm; NFe3O4s -= NFe3O4m; NSiO2s -= NSiO2m;
            MMatte =MAppMatte + MSlagInMatte ;
            MSlag = MAppSlag + MMatteInSlag;
            double NMatte = NCu2Sm + NCum + NFeSm + NFeOm + NFe3O4m + NAsm + NSbm + NBim + NPbm + NPbSm + NZnSm;
            double NSlag = NCu2Os + NFeOs + NFe3O4s + NCu2Ss + NCus + NFeSs + NAs2O3s + NSb2O3s + NBi2O3s + NPbOs + NZnOs + NCaOs + NMgOs + NAl2O3s + NSiO2s;
            double NGas = NO2g + NSO3g + NS2g + NSO2g + NAs2g + NAsSg + NAsOg + NSbOg + NSbSg + NBi2O3g + NPbOg + NPbSg + NZng + NZnOg + NN2g + NH2Og;
            
            double CuMatte = (63.5 * 2 *NCu2Sm +63.5*NCum) / MMatte* 100;
            double CuSlag = (63.5 * 2 * (NCu2Ss + NCu2Os) + 63.5 * NCus) / MSlag * 100;
            double FeOMatte = (56 + 16) * NFeOm / MMatte*100;
            double Fe3O4Matte = (56 * 3 + 16 * 4) * NFe3O4m / MMatte*100;
            double FeOSlag = (56 + 16) * NFeOs / MSlag *100;
             Fe3O4Slag = (56 * 3 + 16 * 4)*NFe3O4s / MSlag *100;
            double Cu2SSlag = (63.5 * 2 + 32)*NCu2Ss / MSlag*100;
            double Cu2OSlag = (63.5 * 2 + 16) * NCu2Os / MSlag * 100;
            double SiO2Matte = (28 + 2 * 16) * NSiO2m / MMatte * 100;
            double SiO2Slag = (28 + 2 * 16) * NSiO2s / MSlag * 100;
            double FeSiO2 = (56 * NFeOs + 56 * 3 * NFe3O4s) / (NSiO2s * (28 + 2 * 16));
            double SO2Gas = NSO2g / NGas;
            double S2Gas = NS2g / NGas;
            double H2OGas = NH2Og / NGas;
            double FeMatte = (56 * NFeSm + 56 * 3 * NFe3O4m +56 * NFeOm) / MMatte * 100;
            double FeSlag = (56 * NFeSs + 56 * 3 * NFe3O4s + 56 * NFeOs) / MSlag * 100;
            double SMatte = (32 * NCu2Sm + 32 * NFeSm + 32 * NPbSm + 32 * NZnSm) / MMatte * 100;
            double SSlag = (32 * NCu2Ss + 32 * NFeSs) / MSlag*100;
            //As
            double AsMatte = 100 * NAsm / (NAsm + 2 * NAs2O3s + +NAs2g + NAsSg + NAsOg); 
            double AsSlag =  100 * 2 * NAs2O3s / (NAsm + 2 * NAs2O3s + +NAs2g + NAsSg + NAsOg);
            double AsGas = 100 - AsMatte - AsSlag;
            double AsInMatte = (74.9 * NAsm) / MMatte*100;
            double AsInSlag = (74.9 * 2 * NAs2O3s) / MSlag*100;
            //Sb
            double SbMatte = 100 * NSbm / (NSbm + 2 * NSb2O3s + NSbOg + NSbSg);
            double SbSlag =  100 * 2 * NSb2O3s / (NSbm + 2 * NSb2O3s + NSbOg + NSbSg);
            double SbGas = 100 - SbMatte - SbSlag;
            double SbInMatte = (121.7 * NSbm) / MMatte * 100;
            double SbInSlag = (121.7 * 2 * NSb2O3s) / MSlag*100;
            //Bi
            double BiMatte = 100 * NBim / (NBim + 2 * NBi2O3s + 2*NBi2O3g);
            double BiSlag =  100 * 2 * NBi2O3s / (NBim + 2 * NBi2O3s + 2*NBi2O3g);
            double BiGas = 100 - BiMatte - BiSlag;
            double BiInMatte = (208.9 * NBim) / MMatte * 100;
            double BiInSlag = (208.9 * 2 * NBi2O3s) / MSlag * 100;
            //Pb
            double PbMatte = 100 * (NPbm + NPbSm) / (NPbm + NPbSm + NPbOs + NPbOg + NPbSg);
            double PbSlag = 100 * NPbOs / (NPbm + NPbSm + NPbOs + NPbOg + NPbSg);
            double PbGas = 100 - PbMatte - PbSlag;
            double PbInMatte = (207 * NPbm +207* NPbSm) / MMatte * 100;
            double PbInSlag = (207 * NPbOs) / MSlag * 100;
            //Zn
            double ZnMatte = 100 * NZnSm / (NZnSm + NZnOs + NZng + NZnOg);
            double ZnSlag =  100 * NZnOs / (NZnSm + NZnOs + NZng + NZnOg);
            double ZnGas = 100 - ZnMatte - ZnSlag;
            double ZnInMatte = (65.5 * NZnSm) / MMatte * 100;
            double ZnInSlag = (65.5 * NZnOs) / MSlag * 100;

            
            double[,] Mid = new double[8,4];
            Mid[0, 0] = CuMatte; Mid[0, 1] = CuSlag; Mid[0, 2] = Cu2SSlag; Mid[0, 3] = Cu2OSlag;
            Mid[1, 0] = FeOMatte; Mid[1, 1] = Fe3O4Matte; Mid[1, 2] = FeOSlag; Mid[1, 3] = Fe3O4Slag;
            Mid[2, 0] = FeSiO2; Mid[2, 1] =SiO2Matte ; Mid[2, 2] = SO2Gas; Mid[2, 3] = H2OGas;
            Mid[3, 0] = AsMatte; Mid[3, 1] = AsSlag ; Mid[3, 2] = AsGas;
            Mid[4, 0] = SbMatte; Mid[4, 1] = SbSlag ; Mid[4, 2] = SbGas;
            Mid[5, 0] = BiMatte; Mid[5, 1] = BiSlag ; Mid[5, 2] = BiGas;
            Mid[6, 0] = PbMatte; Mid[6, 1] = PbSlag ; Mid[6, 2] = PbGas;
            Mid[7, 0] = ZnMatte; Mid[7, 1] = ZnSlag; Mid[7, 2] = ZnGas;

            float[,] tempMid = new float[Mid.GetLength(0), Mid.GetLength(1)];
            for (int i = 0; i < Mid.GetLength(0); i++)
                for (int j = 0; j < Mid.GetLength(1); j++)
                {
                    tempMid[i, j] = (float)Mid[i, j];
                }
            return tempMid;
        }

        //Initialization of velocity to make the velocity satisfy Av = 0;
        public Matrix InitVelocity(Matrix A)
        {
            int Row = A.Row;
            int Col = A.Col;
            Matrix v1 = new Matrix(Col - Row); 
            Matrix A1 = new Matrix(Row,Row);  //Seperate matrix A into A1,A2;
            Matrix A2 = new Matrix(Row,Col-Row);
            
            //Get the matrix by the Guass elimilation and get the index of independent colum vector in the Gauss Matrix;
            Matrix Gauss = A.Gauss();
            int[] Index = Gauss.IndexofIndependentMatrix();      //the member in the Index is 0 and 1; 1 represents the independent colum vector in the Gauss; 
            int one = Index.Count(s => s == 1);   //in this application ,Row equals one, Col substract Row that equals zero;
            int zero = Index.Count(s =>s == 0);
            int[] IndexofOne = new int[one];    //record the index of number 1 in the Index;
            int[] IndexofZero = new int[zero];   //record the index of number 0 in the Index;
            int TempOne = 0,TempZero = 0;
            for (int i = 0; i < Index.Length; i++)
            {   
                
                if(Index[i] == 1)
                {
                    IndexofOne[TempOne] = i;
                    TempOne++;
                }
                else
                {
                    IndexofZero[TempZero] = i;
                    TempZero++;
                }
            }
            //set the members in A1 and A2 by A;
            A = A.Transpose();   //due to the next A1 and A2 is assigned by the row of A;
            for (int i = 0; i < one; i++)
            {
                A1.SetMatrixRow(i, 0, Row - 1, IndexofOne[i], 0, Row - 1, A);
            }
            A2 = A2.Transpose();
            for (int i = 0; i < zero; i++)
            {
                A2.SetMatrixRow(i, 0, Row - 1, IndexofZero[i], 0, Row - 1, A);
            }
            A2 = A2.Transpose(); //A = A.Transpose();

            //random assign the v1
            for (int i = 0; i < zero; i++) 
            {
                v1.SetMatrixValue(i,0.5*rd.NextDouble()); 
            }
            Matrix v2 = -A1.InverseMatrix() * A2 * v1.Transpose();
            v2 = v2.Transpose();
            //combine v1 and v2;
            Matrix v = new Matrix(Index.Length);
            TempOne = 0; TempZero =0; // afresh the TempOne and TempZero for v setting;
            for (int i = 0; i < Index.Length; i++)
            { 
                if(Index[i] == 1)
                {
                    v.SetMatrixValue(i, v2.GetMatrixValue(TempOne));
                    TempOne++;
                }
                else
                {
                    v.SetMatrixValue(i, v1.GetMatrixValue(TempZero));
                    TempZero++;
                }
            }
                return v;
        }

        //The function that get the coefficient matrix of substance;
        public Matrix GetA(ArrayList StrList, string Str)
        {
            Matrix A = new Matrix(StrList.Count);
            int z = 0;
            foreach (string List in StrList)
            {

                int temp = List.IndexOf(Str);
                if (char.IsNumber(List, temp + Str.Length))
                {
                    A[0, z++] = int.Parse(List[temp + Str.Length].ToString());
                }
                else
                {
                    A[0, z++] = 1;
                }
            }
            return A;
        }
        public Matrix GetASO(ArrayList StrList)
        {
            Matrix ASO = new Matrix(2, StrList.Count);
            int z = -1, s = -1;
            foreach (string List in StrList)
            {
                z++; s++;
                for (int i = 0; i < List.Length - 1; i++)
                {
                    if (List[i] == 'S' && List[i + 1] != 'b')
                    {
                        if (char.IsNumber(List[i + 1]))
                        {
                            ASO[0, z] = int.Parse(List[i + 1].ToString());
                        }
                        else { ASO[0, z] = 1; }
                    }
                }

                int temp = List.IndexOf("O");
                if (temp > -1)
                {
                    if (char.IsNumber(List, temp + 1))
                    {
                        ASO[1, s] = int.Parse(List[temp + 1].ToString());
                    }
                    else
                    {
                        ASO[1, s] = 1;
                    }
                }
            }
            return ASO;
        }

        //Used to compare with the double 0;
        private double esp1 = Math.Pow(10, -10);

        //Used for the control of precision and decide the updating step;
        private double esp2 = Math.Pow(10, -30);

       private InputInfo InputPara;
       public InputInfo INPUTPARA
        {
            get { return InputPara; }
            set { InputPara = value; }
        }

        
        public delegate double Myfunc(Matrix x,InputInfo InputPara);
        //Random seed;
        public Random rd = new Random();

        //Used as middle parameter to deliver to the GibbsFunc function;
        private double[] MidVariable;
        public double[] MIDVARIABLE
        {
            get { return MidVariable; }
            set { MidVariable = value; }
        }

        private double[] AdaptiveValue;
        public double[] ADAPTIVEVALUE
        {
            get { return AdaptiveValue; }
            set { AdaptiveValue = value; }
        }

        private double[,] AdaptiveResults;
        public double[,] ADAPTIVERESULTS
        {
            get { return AdaptiveResults; }
            set { AdaptiveResults = value; }
        }
     
        private double[] ConstrainErr;
        public double[] CONSTRAINERR
        {
            get { return ConstrainErr; }
            set { ConstrainErr = value; }
        }
    }
}
