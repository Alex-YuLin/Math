using System;
using System.Collections.Generic;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;
using Math.define;
using Math.Error;

namespace Math
{
    public class Math
    {
        ///////////////////////////////
        ///         Para            ///
        ///////////////////////////////


        ///////////////////////////////
        ///         Sub             ///
        ///////////////////////////////


        private double To360(double theda)
        {
            double O_theda = 0;
            if (theda < 0)
                //O_theda = (180 + (180 - theda))>360? 180 + (180 + theda) : 180 + (180 - theda) ;
                O_theda = 180 + (180 + theda);
            else
                O_theda = theda;

            return O_theda;
        }

        private double RadToDeg(double radi)
        {
            return radi * (180 / System.Math.PI);
        }

        private double DegOfP1toP2(Dot p1, Dot p2)
        {
            double Deg = 0;
            double dx = p2.pX - p1.pX;
            double dy = p2.pY - p1.pY;
            Deg = To360(RadToDeg(System.Math.Atan2(dy, dx)));

            Log.Pushlist_debug("deg: " + RadToDeg(System.Math.Atan2(dy, dx)));
            Log.Pushlist_debug("To360: " + To360(RadToDeg(System.Math.Atan2(dy, dx))));

            return Deg;
        }

        private double LineLenth(Dot p1, Dot p2)
            {
                double pLenth = 0;

                double dx = System.Math.Abs(p1.pX - p2.pX);
                double dy = System.Math.Abs(p1.pY - p2.pY);
                pLenth = System.Math.Sqrt(dx * dx + dy * dy);

                return pLenth;
            }

        public double GetValue(double inputW, double[] W, double[] P)
        {
            if (W.Length != P.Length) return -1;
            int itemNum = W.Length;
            #region Myself
            double value = 0.0;
            //如果初始的离散点为空, 返回0
            if (itemNum < 1) { return value; }
            //如果初始的离散点只有1个, 返回该点对应的Y值
            if (itemNum == 1) { value = W[0]; return value; }
            //如果初始的离散点只有2个, 进行线性插值并返回插值

            int W_start = 0, W_end = 0;

            for (int i = 0; i < W.Length; i++)
            {
                if (W[i] <= inputW)
                    W_start = i;
                else if (W[i] > inputW)
                {
                    W_end = i;
                    break;
                }
            }

            double a1 = (inputW - W[W_start]) * (P[W_end] - P[W_start]);
            double a2 = (W[W_end] - W[W_start]);
            double a3 = P[W_start];
            if (a2 == 0) a2 = 1;
            decimal out1 = (decimal)a1 / (decimal)a2 + (decimal)a3;
            value = Convert.ToDouble(out1);



            #endregion

            return value;
        }


        ///////////////////////////////
        ///         Meth            ///
        ///////////////////////////////

        #region Method

        public double CalRotateDegree(Dot CG_Draw, Dot CG_Table, List<Dot> Draw_Points, List<Dot> Table_Points)
            {
                Log.Pushlist_debug("Start Cal RotetaDegree");
                double pDegree = 0;
                try
                {
                    // define
                    List<double> DegreeDiffs = new List<double>();
                    List<double> DefineCWorCCw = new List<double>();

                    // chk
                    if (Table_Points.Count < 1) throw new Exception("Points < 1");

                    // Cal CG and Points degree offset
                    for (int i = 0; i < Table_Points.Count; i++)
                    {
                        double DraDegree = DegOfP1toP2(CG_Draw, Draw_Points[i]);
                        double ActDegree = DegOfP1toP2(CG_Table, Table_Points[i]);

                        Log.Pushlist_debug("Draw Degree = " + DraDegree);
                        Log.Pushlist_debug("Table Degree = " + ActDegree);

                        double CCWDegDiff = (ActDegree >= DraDegree ? ActDegree - DraDegree : (360.0 - DraDegree) + ActDegree);


                        DegreeDiffs.Add(CCWDegDiff);
                    }

                    // AverageAngleDiff
                    uint DegDiffIndex = 0;
                    uint LargestDegDiffIndex = 0;
                    double LastDegDiff = 0.0;
                    double LargestDiffOfDegDiff = 0.0;
                    for (int i = 0; i < DegreeDiffs.Count; i++)
                    {
                        if (i == 0)
                        {
                            LargestDegDiffIndex = 0;
                            LastDegDiff = DegreeDiffs[i];
                            LargestDiffOfDegDiff = 0.0;
                        }
                        else
                        {
                            if ((DegreeDiffs[i] - LastDegDiff) > LargestDiffOfDegDiff)
                            {
                                LargestDegDiffIndex = DegDiffIndex;
                                LargestDiffOfDegDiff = DegreeDiffs[i] - LastDegDiff;
                            }

                            LastDegDiff = DegreeDiffs[i];
                        }

                        ++DegDiffIndex;
                    }

                    //  計算尾到頭(注意全部 degdiff 都是零度的狀況)
                    {
                        double LastDegDiffToAxisX = (LastDegDiff == 0.0 ? 0.0 : (360.0 - LastDegDiff));
                        double LastDiffOfDefDiff = LastDegDiffToAxisX + DegreeDiffs[0];

                        if (LastDiffOfDefDiff > LargestDiffOfDegDiff)
                        {
                            LargestDegDiffIndex = 0;
                            LargestDiffOfDegDiff = LastDiffOfDefDiff;
                        }
                    }

                    Log.Pushlist_debug("鈍角 (360 - 銳角) = " + LargestDiffOfDegDiff.ToString());


                    double DegAdj = 360 - DegreeDiffs[(int)LargestDegDiffIndex];

                    Log.Pushlist_debug("DegDiff 調整值 = " + DegAdj.ToString());

                    // 計算調整後的DegDiff
                    List<double> DegreeDiffSetOutput = new List<double>();
                    for (int i = 0; i < DegreeDiffs.Count; i++)
                    {
                        double D = DegreeDiffs[i] + DegAdj;

                        //D = (D - 360)<0.05 && (D - 360 )> 0 ? D - 360 : D;

                        Log.Pushlist_debug("調整後 DegDiff:" + D.ToString());

                        DegreeDiffSetOutput.Add(D);

                    }

                    // 計算平均值
                    double AverageAngleDiff = 0;
                    foreach (var om in DegreeDiffSetOutput)
                    {
                        AverageAngleDiff += om;
                    }

                    AverageAngleDiff = AverageAngleDiff / DegreeDiffSetOutput.Count;
                    AverageAngleDiff = AverageAngleDiff + (360 - DegAdj);
                    if (AverageAngleDiff >= 360) AverageAngleDiff -= 360;

                    pDegree = AverageAngleDiff;


                    //// 判斷是CW or CCW
                    //foreach (var om in DegreeDiffs)
                    //{
                    //    if (System.Math.Abs(om) > 180) 
                    //    pDegree = om > 0 ? -pDegree : pDegree;
                    //    break;
                    //}


                }
                catch (Exception e)
                {
                    Log.Pushlist(-1, e.Message);
                    return -1;
                }
                return pDegree;
            }

        public double CalScale(Dot CG_Draw, Dot CG_Table, List<Dot> Draw_Points, List<Dot> Table_Points)
        {
            Log.Pushlist_debug("Start Cal Scale");
            double pScale = 0;
            try
            {
                // chk 
                if (Table_Points.Count < 1) return -1;

                //
                double Scale = 0.0; // 後面是用加的

                // 定位點至Working Home 的平均距離
                double AverageDisToActualWorkingHome = 0;
                double AverageDisToDrawingWorkingHome = 0;

                double ActualSum = 0;
                double DrawingSum = 0;

                for (int i = 0; i < Table_Points.Count; i++)
                {
                    AverageDisToActualWorkingHome += LineLenth(CG_Table, Table_Points[i]);
                    AverageDisToDrawingWorkingHome += LineLenth(CG_Draw, Draw_Points[i]);


                    ActualSum += AverageDisToActualWorkingHome;
                    DrawingSum += AverageDisToDrawingWorkingHome;

                    double ScaleThisPoint = 1.0;

                    if (DrawingSum == 0)
                    {
                        Log.Pushlist_debug("DrawingSum == 0, 此單一次計算比例設為 1 ");

                        ScaleThisPoint = 1.0;
                    }
                    else
                    {
                        ScaleThisPoint = ActualSum / DrawingSum;
                    }

                    Log.Pushlist_debug("本次計算Scale = " + ScaleThisPoint.ToString());

                    Scale += ScaleThisPoint;

                }

                Scale /= Table_Points.Count;

                Log.Pushlist_debug("全圖縮放比例 Scale = " + Scale);

                pScale = Scale;



            }
            catch (Exception e)
            {
                Log.Pushlist(-2, e.Message);
                return -1;
            }
            return pScale;

        }

        public Dot CalScaleXY(double pDegree, Dot CG_Draw, Dot CG_Table, List<Dot> Draw_Points, List<Dot> Table_Points)
        {
            Log.Pushlist_debug("Start Cal ScaleXY");
            Dot pScale = new Dot() { pX = 0, pY = 0 };
            try
            {
                // chk 
                if (Table_Points.Count < 1) return new Dot() { pX = 0, pY = 0 };

                //
                double RotateDegree = pDegree;

                double ScaleX = 0.0; // 後面是用加的
                double ScaleY = 0.0; // 後面是用加的

                Dot RotatedActualPos = new Dot() { pX = 0, pY = 0 };

                double ActualDisToHomeX = 0;
                double ActualDisToHomeY = 0;
                double DrawingDisToHomeX = 0;
                double DrawingDisToHomeY = 0;
                double ScaleXThisPoint = 1;
                double ScaleYThisPoint = 1;

                for (int i = 0; i < Table_Points.Count; i++)
                {
                    RotatedActualPos.pX = Table_Points[i].pX;
                    RotatedActualPos.pY = Table_Points[i].pY;

                    // 轉回去
                    RotatedActualPos = Rotate(CG_Table, RotatedActualPos, RotateDegree);

                    // XY 分量差
                    ActualDisToHomeX = System.Math.Abs(RotatedActualPos.pX - CG_Table.pX);
                    ActualDisToHomeY = System.Math.Abs(RotatedActualPos.pY - CG_Table.pY);
                    DrawingDisToHomeX = System.Math.Abs(Draw_Points[i].pX - CG_Draw.pX);
                    DrawingDisToHomeY = System.Math.Abs(Draw_Points[i].pY - CG_Draw.pY);

                    // 計算單次Scale
                    if (DrawingDisToHomeX <= 5)
                    {
                        Log.Pushlist_debug("log : DrawingDisToHomeX <=5 ,設定此分量 ScaleX = 1.0 ");

                        ScaleXThisPoint = 1;
                    }
                    else
                    {
                        ScaleXThisPoint = ActualDisToHomeX / DrawingDisToHomeX;
                    }
                    if (DrawingDisToHomeY <= 5)
                    {
                        Log.Pushlist_debug("log : DrawingDisToHomeX <=5 ,設定此分量 ScaleY = 1.0");

                        ScaleYThisPoint = 1;
                    }
                    else
                    {
                        ScaleYThisPoint = ActualDisToHomeY / DrawingDisToHomeY;
                    }


                    Log.Pushlist_debug("此次 ScaleX = " + ScaleXThisPoint.ToString());
                    Log.Pushlist_debug("此次 ScaleY = " + ScaleYThisPoint.ToString());

                    ScaleX += ScaleXThisPoint;
                    ScaleY += ScaleYThisPoint;

                }

                ScaleX /= Table_Points.Count;
                ScaleY /= Table_Points.Count;


                Log.Pushlist_debug("全圖 ScaleX = " + ScaleX.ToString());
                Log.Pushlist_debug("全圖 ScaleY = " + ScaleY.ToString());

                pScale = new Dot() { pX = ScaleX, pY = ScaleY };
            }
            catch (Exception e)
            {
                Log.Pushlist(-2, e.Message);
                return new Dot() { pX = 0, pY = 0 };

            }
            return pScale;

        }

        public Dot CalCG(List<Dot> db)
        {
            Dot tamp = new Dot();
            foreach (var om in db)
            {
                tamp.pX += om.pX;
                tamp.pY += om.pY;
            }
            tamp.pX /= db.Count;
            tamp.pY /= db.Count;
            return tamp;
        }

        /// <summary>
        /// anticlockwise
        /// </summary>
        /// <param name="center"></param>
        /// <param name="p1"></param>
        /// <param name="degree"></param>
        /// <returns></returns>
        public Dot Rotate(Dot center, Dot p1, double degree)
        {
            Dot tmp = new Dot();
            double angleHude = degree * System.Math.PI / 180;/*角度變成弧度*/
            double x1 = (p1.pX - center.pX) * System.Math.Cos(angleHude) + (p1.pY - center.pY) * System.Math.Sin(angleHude) + center.pX;
            double y1 = -(p1.pX - center.pX) * System.Math.Sin(angleHude) + (p1.pY - center.pY) * System.Math.Cos(angleHude) + center.pY;
            tmp.pX = (float)x1;
            tmp.pY = (float)y1;
            return tmp;
        }

        public bool GetCross(PointF lineFirstStar, PointF lineFirstEnd, PointF lineSecondStar, PointF lineSecondEnd, out bool bolFind, out PointF PIntersection)
            {
                /*
                 * L1，L2都存在斜率的情況：
                 * 直線方程L1: ( y - y1 ) / ( y2 - y1 ) = ( x - x1 ) / ( x2 - x1 ) 
                 * => y = [ ( y2 - y1 ) / ( x2 - x1 ) ]( x - x1 ) + y1
                 * 令 a = ( y2 - y1 ) / ( x2 - x1 )
                 * 有 y = a * x - a * x1 + y1   .........1
                 * 直線方程L2: ( y - y3 ) / ( y4 - y3 ) = ( x - x3 ) / ( x4 - x3 )
                 * 令 b = ( y4 - y3 ) / ( x4 - x3 )
                 * 有 y = b * x - b * x3 + y3 ..........2
                 * 
                 * 如果 a = b，則兩直線平等，否則， 聯解方程 1,2，得:
                 * x = ( a * x1 - b * x3 - y1 + y3 ) / ( a - b )
                 * y = a * x - a * x1 + y1
                 * 
                 * L1存在斜率, L2平行Y軸的情況：
                 * x = x3
                 * y = a * x3 - a * x1 + y1
                 * 
                 * L1 平行Y軸，L2存在斜率的情況：
                 * x = x1
                 * y = b * x - b * x3 + y3
                 * 
                 * L1與L2都平行Y軸的情況：
                 * 如果 x1 = x3，那麼L1與L2重合，否則平等
                 * 
                */
                bool retrr = true;
                //Point PTemp = new Point(0, 0);
                float a = 0, b = 0;
                int state = 0;
                bolFind = false;
                PointF PInterSection = new PointF();

                if (lineFirstStar.X != lineFirstEnd.X)
                {
                    a = (lineFirstEnd.Y - lineFirstStar.Y) / (lineFirstEnd.X - lineFirstStar.X);
                    state |= 1;
                }
                if (lineSecondStar.X != lineSecondEnd.X)
                {
                    b = (lineSecondEnd.Y - lineSecondStar.Y) / (lineSecondEnd.X - lineSecondStar.X);
                    state |= 2;
                }
                switch (state)
                {
                    case 0: //L1與L2都平行Y軸
                        {
                            if (lineFirstStar.X == lineSecondStar.X)
                            {
                                //throw new Exception("兩條直線互相重合，且平行於Y軸，無法計算交點。");                    
                                bolFind = false;
                            }
                            else
                            {
                                //throw new Exception("兩條直線互相平行，且平行於Y軸，無法計算交點。");                      
                                bolFind = false;
                            }
                        }
                        break;
                    case 1: //L1存在斜率, L2平行Y軸
                        {
                            float x = lineSecondStar.X;
                            float y = (lineFirstStar.X - x) * (-a) + lineFirstStar.Y;
                            bolFind = false;
                            break;
                        }
                    case 2: //L1 平行Y軸，L2存在斜率
                        {
                            float x = lineFirstStar.X;
                            //網上有相似代碼的，這一處是錯誤的。你可以對比case 1 的邏輯 進行分析
                            //源code:lineSecondStar * x + lineSecondStar * lineSecondStar.X + p3.Y;
                            float y = (lineSecondStar.X - x) * (-b) + lineSecondStar.Y;
                            PInterSection.X = (int)x;
                            PInterSection.Y = (int)y;
                            bolFind = true;
                            break;
                        }
                    case 3: //L1，L2都存在斜率
                        {
                            if (a == b)
                            {
                                // throw new Exception("兩條直線平行或重合，無法計算交點。");
                                bolFind = false;
                            }
                            else
                            {
                                float x = (a * lineFirstStar.X - b * lineSecondStar.X - lineFirstStar.Y + lineSecondStar.Y) / (a - b);
                                float y = a * x - a * lineFirstStar.X + lineFirstStar.Y;
                                //return new PointF(x, y);
                                PInterSection.X = (int)x;
                                PInterSection.Y = (int)y;
                                bolFind = true;
                            }
                            break;
                        }
                }
                // throw new Exception("不可能發生的情況");
                PIntersection = PInterSection;//out 
                                              //PInterSection = PTemp;
                return retrr;
            }

        public bool PowerConvert(string _path, double W, out double value)
        {
            double tamp = 0;
            value = tamp;
            try
            {
                #region 讀取雷射Power Table

                string path = _path;
                StreamReader strr2 = new StreamReader(path);
                string ori2 = strr2.ReadToEnd();
                strr2.Close();
                string[] jj = ori2.Split(new string[1] { "\r\n" }, StringSplitOptions.RemoveEmptyEntries);
                double[] Laser_W = new double[jj.Length];
                double[] Laser_Per = new double[jj.Length];
                for (int i = 0; i < jj.Length; i++)
                {
                    string[] kk = jj[i].Split(',');
                    Laser_Per[i] = Convert.ToDouble(kk[1]);
                    Laser_W[i] = Convert.ToDouble(kk[0]);
                }

                #endregion

                #region 換算出來的值

                tamp = GetValue(W, Laser_W, Laser_Per);

                #endregion
            }
            catch (Exception ee)
            {
                return false;
            }
            value = tamp;
            return true;
        }

        #endregion
    }
}

namespace Math.define
{
    public struct Dot
    {
        public Dot(double x, double y)
        {
            pX = x;
            pY = y;
        }
        public double pX;
        public double pY;
    }
}

namespace Math.Distortion.Defind
{

}
namespace Math.Error
{
    class Log
    {
        public static void Pushlist(int NGCode, string what
                                    , [System.Runtime.CompilerServices.CallerMemberName] string Who = ""
                                    , [CallerLineNumber] int line = 0
                                    , [CallerFilePath] string path = "")
        {
            //初始化字串
            string NGString = "";
            string NGtime = DateTime.Now.ToString("yyyy/MM/dd HH:mm:ss");


            string[] ori = path.Split(new string[1] { "\\" }, StringSplitOptions.RemoveEmptyEntries);
            string FileN = ori[ori.Length - 1];
            // 重組字串
            NGString = "# " + NGtime + ", " + FileN + " > func:" + Who + ", " + "line: " + line + ", " + NGCode.ToString() + ", " + what;
            Console.WriteLine(NGString);

            ////寫入文檔
            //string pat = Application.StartupPath + "\\Log\\Math_Distortion.log";
            //StreamWriter strr = new StreamWriter(pat, true);
            //strr.WriteLine(NGString, false);
            //strr.Close();

        }


        public static void Pushlist_debug(string what
                                            , [CallerMemberName] string who = ""
                                            , [CallerLineNumber] int line = 0
                                            , [CallerFilePath] string path = "")
        {
#if !debug
            return;
#endif
            //初始化字串
            string NGString = "";
            string NGtime = DateTime.Now.ToString("yyyy/MM/dd HH:mm:ss");

            string[] ori = path.Split(new string[1] { "\\" }, StringSplitOptions.RemoveEmptyEntries);
            string FileN = ori[ori.Length - 1];
            // 重組字串
            NGString = "# " + NGtime + ", " + FileN + " > func:" + who + ", " + "line: " + line + ", " + what;

            Console.WriteLine(NGString);

#if debug
        //寫入文檔
        //string pat = Application.StartupPath + "\\Log\\Math_Distortion.log";
        //StreamWriter strr = new StreamWriter(pat, true);
        //strr.WriteLine(NGString, false);
        //strr.Close();
#endif
        }
    }
}


