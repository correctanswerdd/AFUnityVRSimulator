using System.Collections;
using System;
using System.Collections.Generic;
using UnityEngine;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using System.Linq;
using Random=UnityEngine.Random;
 
public class MeshBuilder : MonoBehaviour {

    // public Vector3 point = Vector3.up;
    // public int numberOfPoints = 10;

    int[][] LoadFaces(){
        TextAsset txt = Resources.Load("faces") as TextAsset;
        string[] str = txt.text.Split('\n');
        int[][] triangles = new int[str.Length][];
        for(int i=0; i<str.Length; i++)
        {
            string[] ss = str[i].Split(',');
            triangles[i] = new int[]{Int32.Parse(ss[0]), Int32.Parse(ss[1]), Int32.Parse(ss[2])};
        }
        //Debug.Log("Load faces "+str.Length);
        return triangles;
    }

    int[][] LoadFaces_flip(){
        TextAsset txt = Resources.Load("faces") as TextAsset;
        string[] str = txt.text.Split('\n');
        int[][] triangles = new int[str.Length][];
        for(int i=0; i<str.Length; i++)
        {
            string[] ss = str[i].Split(',');
            triangles[i] = new int[]{Int32.Parse(ss[0]), Int32.Parse(ss[2]), Int32.Parse(ss[1])};
        }
        //Debug.Log("Load faces flip"+str.Length);
        return triangles;
    }

    Vector3[] LoadVertices(){
        TextAsset txt = Resources.Load("vertices") as TextAsset;
        string[] str = txt.text.Split('\n');
        Vector3[] vertices = new Vector3[str.Length];
        for(int i=0; i<str.Length; i++)
        {
            string[] ss = str[i].Split(',');
            vertices[i] = new Vector3(float.Parse(ss[0]), float.Parse(ss[1]), float.Parse(ss[2]));
        }
        //Debug.Log("Load vertices "+str.Length);
        return vertices;
    }

    Vector3[] LoadVertices_flip(){
        TextAsset txt = Resources.Load("vertices") as TextAsset;
        string[] str = txt.text.Split('\n');
        Vector3[] vertices = new Vector3[str.Length];
        for(int i=0; i<str.Length; i++)
        {
            string[] ss = str[i].Split(',');
            vertices[i] = new Vector3(float.Parse(ss[0])+100, float.Parse(ss[1])+100, float.Parse(ss[2])+100);
        }
        //Debug.Log("Load vertices flip"+str.Length);
        return vertices;
    }


    int[] LoadColors(int idx, int[][] faces){
        TextAsset txt = Resources.Load("color"+idx.ToString()+".txt") as TextAsset;
        string[] str = txt.text.Split('\n');
        if (str.Length != 4064){
            throw new Exception("wrong color file " + idx);
        }
        //Debug.Log("Load Colors " + str.Length);
        int[] verticesColor = new int[str.Length];
        for(int i=0; i<str.Length; i++){
            verticesColor[i] = int.Parse(str[i]);
        }

        int[] colorList = new int[faces.GetLength(0)];
        for(int i=0; i<faces.GetLength(0); i++){
            colorList[i] = verticesColor[faces[i][0]];
            // colorList[i] = Random.Range(0, 64);
        }

        // int[] colorList = new int[7969];
        // for(int i=0; i<7969; i++){
        //     colorList[i] = Random.Range(0, 64);
        // }

        return colorList;
    }

    int[] LoadColors_v2(int[] verticesColor, int[][] faces)
    {
        int[] colorList = new int[faces.GetLength(0)];
        for (int i = 0; i < faces.GetLength(0); i++)
        {
            colorList[i] = verticesColor[faces[i][0]];
            // colorList[i] = Random.Range(0, 64);
        }

        // int[] colorList = new int[7969];
        // for(int i=0; i<7969; i++){
        //     colorList[i] = Random.Range(0, 64);
        // }

        return colorList;
    }

    public Color RandomColor()
    {
        float r = Random.Range(0f,1f);
        float g = Random.Range(0f,1f);
        float b = Random.Range(0f,1f);
        Color color = new Color(r,g,b);
        return color;
    }

    public Color GenerateColor(float a, float b, float c)
    {
        Color color = new Color(a, b, c);
        return color;
    }

    Material[] SetColors(GameObject temp){
        // Material[] materials = new Material[num_color];
        Material[] materials = temp.transform.parent.GetComponent<MeshRenderer>().materials;
        TextAsset txt = Resources.Load("colormap") as TextAsset;
        string[] str = txt.text.Split('\n');
        for (int i = 0; i < str.Length; i++)
        {
            string[] ss = str[i].Split(',');
            materials[i].color = GenerateColor(float.Parse(ss[0]), float.Parse(ss[1]), float.Parse(ss[2]));
        }
        return materials;
    }

    class LapValues
    {
        private double[][] r_edge;
        private double[][] r_lap;
        public LapValues(int nvertex, double[][] edge, double[][] lap)
        {
            r_edge = new double[nvertex][];
            r_lap = new double[nvertex][];
            for (int i = 0; i < nvertex; i++)
            {
                r_edge[i] = new double[nvertex];
                for (int j = 0; j < nvertex; j++)
                {
                    r_edge[i][j] = edge[i][j];
                }
            }
            for (int i = 0; i < nvertex; i++)
            {
                r_lap[i] = new double[nvertex];
                for (int j = 0; j < nvertex; j++)
                {
                    r_lap[i][j] = lap[i][j];
                }
            }
        }
        public double[][] getEdge()
        {
            return r_edge;
        }
        public double[][] getLap()
        {
            return r_lap;
        }
    }

    private static List<int> findNonZeroIndexVector(double[] vector)
    {
        List<int> list = new List<int>();
        for (int i = 0; i < vector.Length; i++)
        {
            if (vector[i] != 0)
            {
                list.Add(i);
            }
        }
        return list;
    }

    private static List<int> findPositiveIndexVector(double[] vector)
    {
        List<int> list = new List<int>();
        for (int i = 0; i < vector.Length; i++)
        {
            if (vector[i] > 0)
            {
                list.Add(i);
            }
        }
        return list;
    }

    private static double[] getSubVector(double[] vector, int[] index)
    {
        double[] result = new double[index.Length]; 
        for (int i = 0; i < index.Length; i++)
        {
            result[i] = vector[index[i]];
        }
        return result;
    }

    private static void setSubVector<T>(ref T[] vector, int[] index, T[] value)
    {
        if (index.Length != value.Length)
        {
            throw new Exception("Length of index and value should be the same!");
        }
        for (int i = 0; i < index.Length; i++)
        {
            vector[index[i]] = value[i];
        }
    }

    private static Vector<double> roundVector(Vector<double> vector)
    {
        double[] arr_src = vector.ToArray();
        double[] arr_dst = new double[vector.Count];
        for (int i = 0; i< vector.Count; i++)
        {
            arr_dst[i] = Math.Round(arr_src[i], MidpointRounding.AwayFromZero);
        }
        Vector<double> result  = Vector<double>.Build.DenseOfArray(arr_dst); ;
        return result;
    }


    

    public class AFTriggerDetect
    {
        private int[][] t;  // edges
        private Vector3[] X_data;  // vertices
        private double[][] edges;
        private double[][] lap1;
        private int[] tr;
        private double a;
        private double b;
        private double c1;
        private double c2;
        private double d;
        private double delta;
        private double deltaT;
        private int trid;  // 0~12
        private int n;  //nvertices
        private int T;
        private int lent;
        private double st;
        private int N;  // N = T * lent = 4800 * 1

        // FHN model variables
        private static Vector<double> u;
        private static Vector<double> v;
        private static Vector<double> F;
        private static Vector<double> G;
        private static double[][] lap;

        // steps
        // private int step;

        public AFTriggerDetect(Vector3[] vertex, int[][] faces, int tridd)
        {
            double[] para = new double[6] { 0.12, 0.013, 0.26, 0.1, 1.0, 0.6 };
            a = para[0];
            b = para[1];
            c1 = para[2];
            c2 = para[3];
            d = para[4];
            delta = para[5];

            tr = new int[13]; 
            int[] trr = new int[13] {2410, 2411, 2417, 2418, 2434, 2439, 2444, 2375, 2325, 2241, 2502, 2527, 2575 };
            for (int i = 0; i < 13; i++)
            {
                tr[i] = trr[i] - 1;
            }

            t = faces;
            X_data = vertex;

            LapValues lapValues = mesh_laplacian(vertex, faces);
            edges = lapValues.getEdge();
            lap1 = lapValues.getLap();

            trid = tridd;

            n = vertex.Length;

            T = 4800;
            lent = 1;
            N = T * lent;
            st = T / 3;
            deltaT = 0.5;
        }
        private LapValues mesh_laplacian(Vector3[] vertex, int[][] faces)
        {
            int nvertex = vertex.Length;
            int nfaces = faces.GetLength(0);

            double[][] edge = new double[nvertex][];
            for (int i = 0; i < nvertex; i++)
            {
                edge[i] = new double[nvertex];
            }

            for (int i = 0; i < nfaces; i++)
            {
                Matrix A = DenseMatrix.OfArray(new double[,] {
                { vertex[faces[i][0]][0], vertex[faces[i][0]][1], vertex[faces[i][0]][2]},
                { vertex[faces[i][1]][0], vertex[faces[i][1]][1], vertex[faces[i][1]][2]},
                { vertex[faces[i][2]][0], vertex[faces[i][2]][1], vertex[faces[i][2]][2]}});

                Matrix B = DenseMatrix.OfArray(new double[,] {
                { vertex[faces[i][1]][0], vertex[faces[i][1]][1], vertex[faces[i][1]][2]},
                { vertex[faces[i][2]][0], vertex[faces[i][2]][1], vertex[faces[i][2]][2]},
                { vertex[faces[i][0]][0], vertex[faces[i][0]][1], vertex[faces[i][0]][2]}});

                var Diff = A.Subtract(B);
                var Norm = Diff.PointwisePower(2).RowSums().PointwiseSqrt();

                edge[faces[i][0]][faces[i][1]] = Norm[0];
                edge[faces[i][1]][faces[i][2]] = Norm[1];
                edge[faces[i][2]][faces[i][0]] = Norm[2];

                edge[faces[i][1]][faces[i][0]] = Norm[0];
                edge[faces[i][2]][faces[i][1]] = Norm[1];
                edge[faces[i][0]][faces[i][2]] = Norm[2];
            }
            
            // Using edge to identify nearest vertices, calculate
            // the Laplacian for an irregular mesh
            double[][] lap = new double[nvertex][];
            for (int i = 0; i < nvertex; i++)
            {
                lap[i] = new double[nvertex];
            }
            for (int i = 0; i < nvertex; i++)
            {
                List<int> k = findNonZeroIndexVector(edge[i]);
                // Debug.Log("k" + k[0]);
                if (k.Count == 0)
                {
                    for (int j = 0; j < nvertex; j++)
                    {
                        lap[i][j] = 0;
                    }
                }
                else
                {
                    int ni = k.Count; // the number of neighbours
                    double hi = getSubVector(edge[i], k.ToArray()).Average(); // the average distance to the neighbours
                    
                    // the average inverse distance to the neighbours
                    double[] a = getSubVector(edge[i], k.ToArray());
                    double[] b = new double[a.Length];
                    for (int j = 0; j < a.Length; j++)
                    {
                        b[j] = 1.0 / a[j];
                    }
                    double invhi = b.Average();
                  
                    // Laplacian of vertex itself
                    lap[i][i] = -(4 / hi) * invhi;

                    // Laplacian of direct neighbours
                    double[] b1 = new double[a.Length];
                    for (int j = 0; j < a.Length; j++)
                    {
                        b1[j] = 1.0 / a[j] * (4 / (hi * ni));
                    }
                    setSubVector<double>(ref lap[i], k.ToArray(), b1);
                }
            }

            // return
            LapValues lapValues = new LapValues(nvertex, edge, lap);
            return lapValues;
        }

        public void initializeFHNModel()
        {
            // Initialize FHN model variables
            double[] uu = new double[n];
            u = Vector<double>.Build.DenseOfArray(uu);

            double[] vv = new double[n];
            v = Vector<double>.Build.DenseOfArray(vv);

            double[] FF = new double[n];
            F = Vector<double>.Build.DenseOfArray(FF);

            double[] GG = new double[n];
            G = Vector<double>.Build.DenseOfArray(GG);

            lap = lap1; // initialize 'lap' with 'lap1'
        }


        public double[] stepColor(int nt)
        {
            Debug.Log("nt" + nt);
            if (nt > N)
            {
                throw new Exception("nt out of bounds!");
            }
            double[] Iex = new double[n];
            if (nt % T < 10)
            {
                // Debug.Log("lap296" + lap[295][250]);
                List<int> idx = findPositiveIndexVector(lap[295]); //stimulate SA node every heart cycle 'T'
                double[] b = new double[idx.Count];
                for (int i = 0; i < idx.Count; i++)
                {
                    b[i] = 5;
                }
                setSubVector<double>(ref Iex, idx.ToArray(), b);
                // Debug.Log("idx" + idx[0] + "-" + idx[1]);

                List<int> idx1 = findPositiveIndexVector(lap[tr[trid]]);  // stimulate AF trigger every 1 / 4 heart cycle 'T/4'
                double[] b1 = new double[idx1.Count];
                for (int i = 0; i < idx1.Count; i++)
                {
                    b1[i] = 3;
                }
                setSubVector<double>(ref Iex, idx1.ToArray(), b1);
            }
            else if (nt >= st && nt <= st + 9)
            {
                List<int> idx1 = findPositiveIndexVector(lap[tr[trid]]);  // stimulate AF trigger every 1 / 4 heart cycle 'T/4'
                double[] b1 = new double[idx1.Count];
                for (int i = 0; i < idx1.Count; i++)
                {
                    b1[i] = 3;
                }
                setSubVector<double>(ref Iex, idx1.ToArray(), b1);
            }
            //Debug.Log("Iex" + Iex[2415]);

            // FHN model
            Vector<double> v_Iex = Vector<double>.Build.DenseOfArray(Iex);
            double[,] lapp = new double[n, n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    lapp[i, j] = lap[i][j];
                }
            }
            Matrix<double> m_lap = DenseMatrix.OfArray(lapp);

            F = c1 * u.PointwiseMultiply(u - a).PointwiseMultiply(1 - u) - c2 * u.PointwiseMultiply(v) + v_Iex; // c1*u.*(u-a).*(1-u)-c2*u.*v + Iex;    
            //Debug.Log("F" + F[2415]);
            G = b * (u - d * v);
            u = u + deltaT * (F + delta * m_lap * u);
            v = v + deltaT * G;
            //Debug.Log("u"+u[2415]);

            // color index
            Vector<double> m = 1 + roundVector(63 * u);
            Vector<double> ones = Vector<double>.Build.Dense(n, 1);
            Vector<double> sixtyfour = Vector<double>.Build.Dense(n, 64);
            m = m.PointwiseMaximum(ones); // m = max(m, 1); 
            m = m.PointwiseMinimum(sixtyfour); // min(m, 64);
            //Debug.Log(m[2415]);

            return (m - 1).ToArray();
        }
    }

    void Start () {
        GameObject temp = new GameObject();
        temp.name = "test";
        //获得父物体组件，材质数组挂在父物体上
        temp.transform.parent = this.transform;
        int[][] triangles = LoadFaces();
        Vector3[] vertices = LoadVertices();

        Mesh mesh = new Mesh();
        mesh.vertices = vertices;
        mesh.subMeshCount = triangles.Length;
        for(int i=0; i<mesh.subMeshCount; i++){
            mesh.SetTriangles( triangles[i] , i );
        }

        // MeshFilter
        MeshFilter meshFilter = temp.GetComponent<MeshFilter>();
        if( meshFilter == null ) meshFilter = temp.AddComponent<MeshFilter>();
        meshFilter.mesh = mesh;

        // MeshRenderer
        MeshRenderer meshRenderer = temp.GetComponent<MeshRenderer>();
        if( meshRenderer == null ) meshRenderer = temp.AddComponent<MeshRenderer>();


        ////////////////////////
        // flip obj
        ////////////////////////
        GameObject temp_flip = new GameObject();
        temp_flip.name = "flip";
        temp_flip.transform.parent = this.transform;
        int[][] triangles_flip = LoadFaces_flip();
        // Vector3[] vertices_flip = LoadVertices_flip();

        Mesh mesh_flip = new Mesh();
        // mesh_flip.vertices = vertices_flip;
        mesh_flip.vertices = vertices;
        mesh_flip.subMeshCount = triangles_flip.Length;
        for(int i=0; i<mesh_flip.subMeshCount; i++){
            mesh_flip.SetTriangles( triangles_flip[i] , i );
        }

        // MeshFilter
        MeshFilter meshFilter_flip = temp_flip.GetComponent<MeshFilter>();
        if( meshFilter_flip == null ) meshFilter_flip = temp_flip.AddComponent<MeshFilter>();
        meshFilter_flip.mesh = mesh_flip;

        // MeshRenderer
        MeshRenderer meshRenderer_flip = temp_flip.GetComponent<MeshRenderer>();
        if( meshRenderer_flip == null ) meshRenderer_flip = temp_flip.AddComponent<MeshRenderer>();


        // Materials
        // Material[] materials_get = tem.transform.parent.GetComponent<MeshRenderer>().materials;
        Material[] materials_get = SetColors(temp);
        Material[] materials_get_flip = SetColors(temp_flip);
        AFTriggerDetect aftrigger = new AFTriggerDetect(vertices, triangles, 0);
        aftrigger.initializeFHNModel();

        StartCoroutine(ChangeColor(triangles, mesh, mesh_flip, meshRenderer, meshRenderer_flip, materials_get, materials_get_flip, aftrigger));
    }

    public IEnumerator ChangeColor(int[][] faces, Mesh mesh, Mesh mesh_flip, MeshRenderer meshRenderer, MeshRenderer meshRenderer_flip, Material[] materials_get, Material[] materials_get_flip, AFTriggerDetect aftrigger)
    {
        var index = 1;
        while (true)
        {
            double[] v_colors = aftrigger.stepColor(index);
            int[] color_idx = v_colors.Select(d => (int)d).ToArray(); 
            ShowGameObject(faces, mesh, meshRenderer, materials_get, color_idx); // convert a double array to integer one
            ShowGameObject(faces, mesh_flip, meshRenderer_flip, materials_get_flip, color_idx);
            index ++;
            if (index >= 4800)
            {
                // index = 1;
                break;
            }

            yield return new WaitForSeconds(0.2f);
        }
    }

    public void ShowGameObject(int[][] faces, Mesh mesh, MeshRenderer meshRenderer, Material[] materials_get, int[] verticesColor)
    {
        int[] color_idx = LoadColors_v2(verticesColor, faces);
        Material[] materials_push = new Material[mesh.subMeshCount];
        for (int i = 0; i < mesh.subMeshCount; i++)
        {
            materials_push[i] = materials_get[color_idx[i]];
        }
        meshRenderer.materials = materials_push;//添加多个材质
    }


    // Update is called once per frame
    void Update () {
     
    }

}