using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MathematicsFactory
{
    public class Delaunay
    {
        public class Point
        {
            public double[] x = new double[3];
            public Point(double x0 = 0, double x1 = 0, double x2 = 0)
            { x[0] = x0; x[1] = x1; x[2] = x2; }
            public Point(Doub x0, Doub x1)
            { x[0] = (double)x0; x[1] = (double)x1; }
            public Point(Point p)
            { for (int i = 0; i < 3; i++) x[i] = p.x[i]; }
            public static bool operator ==(Point a, Point b)
            {
                for (int i = 0; i < 3; i++) if (a.x[i] != b.x[i]) return false;
                return true;
            }
            public static bool operator !=(Point a, Point b)
            {
                for (int i = 0; i < 3; i++) if (a.x[i] != b.x[i]) return true;
                return false;
            }
            public static double Dist(Point a, Point b)
            {
                double dd = 0;
                for (int i = 0; i < 2; i++) dd += SQR(a.x[i] - b.x[i]);
                return Sqrt(dd);
            }
            static double Sqrt(double x)
            { return Math.Sqrt(x); }
            static double SQR(double x)
            { return x * x; }
        }
        public class Triangle
        {
            public List<int> pointlist = new List<int>();
            public List<int> edgelist = new List<int>();
        }
        public class Edge
        {
            public List<int> pointlist = new List<int>();
            public List<int> trilist = new List<int>();
        }
        public class Circle
        {
            public Point Center;
            public double Radius;
        }
        public class TriangleNode
        {
            public int mom; public List<int> daus = new List<int>();
            public bool haschildren = false;
            public Triangle T;
            public TriangleNode(List<int> pointindex)
            {
                T = new Triangle();
                T.pointlist = pointindex.ToList();
            }
        }


        //Triangulation begins
        List<Point> Points;
        List<TriangleNode> Triangles = new List<TriangleNode>();
        List<Edge> Edges = new List<Edge>();
        List<int> newedgeid, newtrisid, newtrisidtemp, edge2legalize, edge2legalizetemp;
        bool PointOnEdge; int edge, ip1, iedge, itri, N, ptid;
        public int[,] TRI, EDGE;
        public Delaunay(Matrix<Doub> XYZSurface)
        {
            double xmin = double.PositiveInfinity, ymin = double.PositiveInfinity, x,
                   xmax = double.NegativeInfinity, ymax = double.NegativeInfinity, y;
            Points = new List<Point>();
            for (int i = 0; i < XYZSurface.Rows; i++)
            {
                x = (double)XYZSurface[i, 0]; y = (double)XYZSurface[i, 1]; Points.Add(new Point(x, y));
                xmin = Math.Min(xmin, x); xmax = Math.Max(xmax, x);
                ymin = Math.Min(ymin, y); ymax = Math.Max(ymax, y);
            }
            Points = Points.Distinct().ToList(); Circle C = new Circle(); N = Points.Count;
            C.Center = new Point(0.5 * (xmin + xmax), 0.5 * (ymin + ymax)); // Creating the center
            Point lo = new Point(xmin, ymin), hi = new Point(xmax, ymax); // Highest and lowest point
            C.Radius = Point.Dist(C.Center, lo); // getting the radius
            double r = C.Radius * 1.001; // Increasing the radius so no point fall on its circumfrence
            //Adding the points of the big triangle to the  list of points. Counter clockwise
            Points.Add(new Point(C.Center.x[0], C.Center.x[1] + 2 * r));
            Points.Add(new Point(C.Center.x[0] - Math.Sqrt(3) * r, C.Center.x[1] - r));
            Points.Add(new Point(C.Center.x[0] + Math.Sqrt(3) * r, C.Center.x[1] - r));

            Triangles.Add(new TriangleNode(new List<int>() { N, N + 1, N + 2 })); Triangles.Last().T.edgelist.AddRange(new int[] { 0, 1, 2 });
            Edges.Add(new Edge()); Edges.Last().pointlist = new List<int>() { N, N + 1 }; Edges.Last().trilist.Add(0);
            Edges.Add(new Edge()); Edges.Last().pointlist = new List<int>() { N + 1, N + 2 }; Edges.Last().trilist.Add(0);
            Edges.Add(new Edge()); Edges.Last().pointlist = new List<int>() { N + 2, N }; Edges.Last().trilist.Add(0);
            iedge = 2; itri = 0;
            for (ptid = 0; ptid < N; ptid++)
            {
                newtrisid = new List<int>(); newedgeid = new List<int>();
                int ParentId = FindContainingTriangle(ptid);
                if (PointOnEdge)
                {
                    List<int> Quadedgelist, Quadpointlist;
                    EdgesAndPoints(edge, out Quadedgelist, out Quadpointlist);
                    Edges.Add(new Edge()); Edges.Last().pointlist = new List<int>() { ptid, Quadpointlist[0] }; newedgeid.Add(++iedge);
                    Edges.Add(new Edge()); Edges.Last().pointlist = new List<int>() { ptid, Quadpointlist[1] }; newedgeid.Add(++iedge);
                    Edges.Add(new Edge()); Edges.Last().pointlist = new List<int>() { ptid, Quadpointlist[2] }; newedgeid.Add(++iedge);
                    Edges.Add(new Edge()); Edges.Last().pointlist = new List<int>() { ptid, Quadpointlist[3] }; newedgeid.Add(++iedge);
                    for (int i = 0; i < 2; i++)
                    {
                        ParentId = Edges[edge].trilist[i];
                        TriangleNode Parent = Triangles[ParentId];
                        for (int j = 0; j < 2; j++)
                        {
                            int m = 2 * i + j, mp1 = (m + 1) % 4;
                            Triangles.Add(new TriangleNode(new List<int>() { ptid, Quadpointlist[m], Quadpointlist[mp1] }));
                            Triangles.Last().mom = ParentId; Parent.daus.Add(++itri); newtrisid.Add(itri);
                            Triangles.Last().T.edgelist = new List<int>() { newedgeid[m], newedgeid[mp1] };
                            foreach (int ii in Quadedgelist)
                            {
                                AddTriangleandEdge(ii, itri);
                                RemoveTriangleandEdge(ii, ParentId);
                            }
                        }
                        Triangles[ParentId].haschildren = true;
                    }
                    edge2legalize = Quadedgelist;
                }
                else
                {
                    TriangleNode Parent = Triangles[ParentId];
                    List<int> Triedgelist = Parent.T.edgelist.ToList(), Tripointlist = Parent.T.pointlist.ToList();
                    Edges.Add(new Edge()); Edges.Last().pointlist = new List<int>() { ptid, Parent.T.pointlist[0] }; newedgeid.Add(++iedge);
                    Edges.Add(new Edge()); Edges.Last().pointlist = new List<int>() { ptid, Parent.T.pointlist[1] }; newedgeid.Add(++iedge);
                    Edges.Add(new Edge()); Edges.Last().pointlist = new List<int>() { ptid, Parent.T.pointlist[2] }; newedgeid.Add(++iedge);
                    for (int i = 0; i < 3; i++)
                    {
                        ip1 = (i + 1) % 3;
                        Triangles.Add(new TriangleNode(new List<int>() { ptid, Tripointlist[i], Tripointlist[ip1] }));
                        Triangles.Last().mom = ParentId; Parent.daus.Add(++itri); newtrisid.Add(itri);
                        Triangles.Last().T.edgelist = new List<int>() { newedgeid[i], newedgeid[ip1] };
                        Edges[newedgeid[i]].trilist.Add(itri); Edges[newedgeid[ip1]].trilist.Add(itri);
                        foreach (int ii in Triedgelist)
                        {
                            AddTriangleandEdge(ii, itri);
                            RemoveTriangleandEdge(ii, ParentId);
                        }
                    }
                    Triangles[ParentId].haschildren = true;
                    edge2legalize = Triedgelist;
                }
                while (edge2legalize.Count > 0)
                {
                    edge2legalizetemp = new List<int>();
                    newtrisidtemp = new List<int>();
                    foreach (int ii in edge2legalize)
                        Legalize(ii);
                    edge2legalize = edge2legalizetemp.ToList();
                    newtrisid = newtrisidtemp.ToList();
                }
                Disp(Triangles.FindAll(xx => !xx.haschildren).ToList());
            }
            List<TriangleNode> Tris = Triangles.FindAll(xx => !xx.haschildren && xx.T.pointlist.All(xxx => xxx < N)).ToList();
            TRI = new int[Tris.Count, 3]; List<int> childrenedges = new List<int>();
            for (int i = 0; i < Tris.Count; i++) for (int j = 0; j < 3; j++)
                { TRI[i, j] = Tris[i].T.pointlist[j]; childrenedges.AddRange(Tris[i].T.edgelist); }
            childrenedges = childrenedges.Distinct().ToList();
            EDGE = new int[childrenedges.Count, 2];
            for (int i = 0; i < childrenedges.Count; i++)
            {
                EDGE[i, 0] = Edges[childrenedges[i]].pointlist[0];
                EDGE[i, 1] = Edges[childrenedges[i]].pointlist[1];
            }
        }

        void Disp(List<TriangleNode> Tris)
        {
            foreach (var Node in Tris)
            {
                foreach (var index in Node.T.pointlist)
                {
                    Console.Write("{0} \t", index);
                }
                Console.WriteLine("");
            }
            Console.WriteLine("");
            Console.WriteLine("");
            Console.WriteLine("");
        }

        

        int FindContainingTriangle(int p)
        {
            int j = 0;
            while (Triangles[j].haschildren)
            {
                for (int i = 0; i < Triangles[j].daus.Count; i++)
                {
                    Triangle T = Triangles[Triangles[j].daus[i]].T;
                    if (IsIn(T, Points[p]))
                    { j = Triangles[j].daus[i]; break; }
                }
            }
            return j;
        }

        void Legalize(int edgeid)
        {
            //Test Legality//
            int tri1 = 0, tri1ind, tri2 = 0; List<int> Quadpointlist = null, Quadedgelist = null, violators = null, flipedgepoint;
            bool IsLegal = false;
            if (Edges[edgeid].pointlist.All(i => i >= N)) // case 1
                IsLegal = true;
            else
            {
                tri1 = Edges[edgeid].trilist[1]; tri2 = Edges[edgeid].trilist[0];
                EdgesAndPoints(edgeid, out Quadedgelist, out Quadpointlist);
                violators = Quadpointlist.FindAll(xx => xx >= N).ToList();
                if (violators.Count == 0 || violators.Count == 2) //case 2 and 4
                {
                    int a = Quadpointlist[0], b = Quadpointlist[1],
                            c = Quadpointlist[2], p = Quadpointlist[3];
                    bool encircle = abcEncirclesp(p, a, b, c);
                    IsLegal = !encircle;
                }
                else  // case 3
                {
                    IsLegal = Edges[edgeid].pointlist.All(xx => xx < N);
                    if (!IsLegal)
                    {
                        double angle = 0; int innerpoint = Edges[edgeid].pointlist.Min();
                        for (int k = 0; k < 2; k++)
                        {
                            TriangleNode Tnode = Triangles[Edges[edgeid].trilist[k]];
                            List<double> angles = ComputeAngles(Edges[edgeid].trilist[k]);
                            int pointindex = Tnode.T.pointlist.IndexOf(innerpoint);
                            angle += angles[pointindex];
                        }
                        IsLegal = angle >= Math.PI;
                    }
                }
            }

            if (!IsLegal)
            {
                //Flipping edge
                tri1 = Edges[edgeid].trilist[1]; tri2 = Edges[edgeid].trilist[0];
                flipedgepoint = new List<int>() { Quadpointlist[1], Quadpointlist[3] };
                Edges.Add(new Edge()); Edges.Last().pointlist = flipedgepoint; iedge++;
                for (int j = 0; j < 2; j++)
                {
                    int jp1 = (j + 1) % 2;
                    Triangles.Add(new TriangleNode(Quadpointlist.FindAll(x => x != Edges[edgeid].pointlist[jp1]).ToList()));
                    Triangles[tri1].daus.Add(++itri); Triangles[tri2].daus.Add(itri); newtrisidtemp.Add(itri);
                    Edges[iedge].trilist.Add(itri); Triangles[itri].T.edgelist.Add(iedge);
                    foreach (int ii in Quadedgelist)
                    {
                        AddTriangleandEdge(ii, itri);
                        RemoveTriangleandEdge(ii, tri1);
                        RemoveTriangleandEdge(ii, tri2);
                    }
                }
                Triangles[tri1].T.edgelist.Remove(edgeid); Triangles[tri2].T.edgelist.Remove(edgeid);
                Triangles[tri1].haschildren = true; Triangles[tri2].haschildren = true;
                edge2legalizetemp.AddRange(Quadedgelist.FindAll(x => !Edges[x].pointlist.Contains(ptid)));
            }
        }

        void EdgesAndPoints(int interfaceedge, out List<int> edges, out List<int> points)
        {
            int i, j, k, l, r = 0, s = 1;
            if (newtrisid.Count > 0) { r = 1; s = 0; }
            points = Triangles[Edges[interfaceedge].trilist[r]].T.pointlist.Union(
                     Triangles[Edges[interfaceedge].trilist[s]].T.pointlist).ToList();
            edges = Triangles[Edges[interfaceedge].trilist[0]].T.edgelist.Union(
                    Triangles[Edges[interfaceedge].trilist[1]].T.edgelist).ToList();
            edges.Remove(interfaceedge);
        }

        bool IsIn(Triangle T, Point P)
        {
            Point p0, p1, p2; double x0, x1, x2, y0, y1, y2;
            Matrix<Doub> A, vp, up; List<Doub> u; bool ans;
            p0 = Points[T.pointlist[0]]; p1 = Points[T.pointlist[1]]; p2 = Points[T.pointlist[2]];
            x0 = p0.x[0]; x1 = p1.x[0]; x2 = p2.x[0]; y0 = p0.x[1]; y1 = p1.x[1]; y2 = p2.x[1];
            A = new Doub[,] { { x1 - x0, x2 - x0 }, { y1 - y0, y2 - y0 } };
            vp = new Doub[,] { { P.x[0] - x0 }, { P.x[1] - y0 } }; up = Matrix<Doub>.Solve(A, vp);
            u = up.ToList(); ans = (u.All(x => x <= 1) && u.All(x => x >= 0) && u.Sum(x => (double)x) <= 1);
            if (ans)
            {
                PointOnEdge = false;
                if (u[0] == 0) { edge = T.edgelist[1]; PointOnEdge = true; }
                if (u[1] == 0) { edge = T.edgelist[0]; PointOnEdge = true; }
                if (u.Sum(x => (double)x) == 1) { edge = T.edgelist[2]; PointOnEdge = true; }
            }
            return ans;
        }

        bool abcEncirclesp(int p, int a, int b, int c)
        {
            Point Pa = Points[a], Pb = Points[b], Pc = Points[c], P = Points[p];
            Circle C = Circumcircle(Pa, Pb, Pc);
            return Point.Dist(C.Center, P) < C.Radius;
        }

        Circle Circumcircle(Point a, Point b, Point c)
        {
            double a0, a1, c0, c1, det, asq, csq, ctr0, ctr1, rad2;
            a0 = a.x[0] - b.x[0]; a1 = a.x[1] - b.x[1];
            c0 = c.x[0] - b.x[0]; c1 = c.x[1] - b.x[1];
            det = a0 * c1 - c0 * a1;
            if (det == 0.0) throw new Exception("no circle three colinear points");
            det = 0.5 / det;
            asq = a0 * a0 + a1 * a1;
            csq = c0 * c0 + c1 * c1;
            ctr0 = det * (asq * c1 - csq * a1);
            ctr1 = det * (csq * a0 - asq * c0);
            rad2 = ctr0 * ctr0 + ctr1 * ctr1;
            Circle ans = new Circle(); ans.Center = new Point(ctr0 + b.x[0], ctr1 + b.x[1]);
            ans.Radius = Math.Sqrt(rad2);
            return ans;
        }

        void AddTriangleandEdge(int iedge, int itri)
        {
            if (Edges[iedge].pointlist.All(xx => Triangles[itri].T.pointlist.Contains(xx)))
            { Triangles[itri].T.edgelist.Add(iedge); Edges[iedge].trilist.Add(itri); }
        }

        void RemoveTriangleandEdge(int iedge, int itri)
        { Triangles[itri].T.edgelist.Remove(iedge); Edges[iedge].trilist.Remove(itri); }

        List<double> ComputeAngles(int triindex)
        {
            List<int> pointlist = Triangles[triindex].T.pointlist;
            List<double> dist = new List<double>(), angs = new List<double>();
            for (int k = 0; k < 3; k++)
            {
                int kp1 = pointlist[(k + 1) % 3], kp2 = pointlist[(k + 2) % 3];
                dist.Add(Point.Dist(Points[kp1], Points[kp2]));
            }
            for (int k = 0; k < 3; k++)
            {
                int kp1 = (k + 1) % 3, kp2 = (k + 2) % 3;
                double a = dist[k], b = dist[kp1], c = dist[kp2],
                       a2 = a * a, b2 = b * b, c2 = c * c, bc = b * c;
                angs.Add(Math.Acos((b2 + c2 - a2) / (2 * bc)));
            }
            return angs;
        }


    }
}
