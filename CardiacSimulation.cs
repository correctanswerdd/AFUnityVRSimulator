using System.Collections;
using System;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;
using UnityEngine.Events;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using System.Linq;
using Random = UnityEngine.Random;

using UnityEngine.XR;
using UnityEngine.XR.Interaction.Toolkit;
using TMPro;

[RequireComponent(typeof(InputData))]

public class CardiacSimulation : MonoBehaviour
{
    int[][] LoadFaces()
    {
        TextAsset txt = Resources.Load("faces") as TextAsset;
        string[] str = txt.text.Split('\n');
        int[][] triangles = new int[str.Length][];
        for (int i = 0; i < str.Length; i++)
        {
            string[] ss = str[i].Split(',');
            triangles[i] = new int[] { Int32.Parse(ss[0]), Int32.Parse(ss[1]), Int32.Parse(ss[2]) };
        }
        Debug.Log("Load faces "+str.Length);
        return triangles;
    }

    int[][] LoadAllFaces()
    {
        TextAsset txt = Resources.Load("faces") as TextAsset;
        string[] str = txt.text.Split('\n');
        int[][] triangles = new int[str.Length * 2][];
        for (int i = 0; i < str.Length; i += 2)
        {
            string[] ss = str[i].Split(',');
            triangles[i] = new int[] { Int32.Parse(ss[0]), Int32.Parse(ss[1]), Int32.Parse(ss[2]) };
            triangles[i + 1] = new int[] { Int32.Parse(ss[0]), Int32.Parse(ss[2]), Int32.Parse(ss[1])};
        }
        Debug.Log("Load all faces " + str.Length * 2);
        return triangles;
    }


    int[][] LoadFaces_flip()
    {
        TextAsset txt = Resources.Load("faces") as TextAsset;
        string[] str = txt.text.Split('\n');
        int[][] triangles = new int[str.Length][];
        for (int i = 0; i < str.Length; i++)
        {
            string[] ss = str[i].Split(',');
            triangles[i] = new int[] { Int32.Parse(ss[0]), Int32.Parse(ss[2]), Int32.Parse(ss[1]) };
        }
        //Debug.Log("Load faces flip"+str.Length);
        return triangles;
    }

    Vector3[] LoadVertices()
    {
        TextAsset txt = Resources.Load("vertices") as TextAsset;
        string[] str = txt.text.Split('\n');
        Vector3[] vertices = new Vector3[str.Length];
        for (int i = 0; i < str.Length; i++)
        {
            string[] ss = str[i].Split(',');
            vertices[i] = new Vector3(float.Parse(ss[0]), float.Parse(ss[1]), float.Parse(ss[2]));
            //Debug.Log("(" + float.Parse(ss[0]) * 0.1f + "," + float.Parse(ss[1]) * 0.1f + "," + float.Parse(ss[2]) * 0.1f + ")");
        }
        Debug.Log("Load vertices " + str.Length);
        return vertices;
    }

    int[] ColorV2F(int[] verticesColor, int[][] faces)
    {
        int[] colorList = new int[faces.GetLength(0)];
        for (int i = 0; i < faces.GetLength(0); i++)
        {
            colorList[i] = verticesColor[faces[i][0]];
        }
        return colorList;
    }

    public Color GenerateColor(float a, float b, float c)
    {
        Color color = new Color(a, b, c);
        return color;
    }

    Material[] SetColors(GameObject temp)
    {
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

    private static int[] findNonZeroIndexVector(double[] vector)
    {
        List<int> list = new List<int>();
        for (int i = 0; i < vector.Length; i++)
        {
            if (vector[i] != 0)
            {
                list.Add(i);
            }
        }
        return list.ToArray();
    }

    private static double avgNeighbourDist(double[] vector, int[] neighbourid, int neighbourN)
    {
        double sum = 0;
        foreach (int i in neighbourid)
            sum += vector[i];
        return sum / neighbourN;
    }

    private static double avgNeighbourDistInverse(double[] vector, int[] neighbourid, int neighbourN)
    {
        double sum = 0;
        foreach (int i in neighbourid)
            sum += 1.0 / vector[i];
        return sum / neighbourN;
    }

    private T[] setSubVector<T>(int[] index, T value)
    {
        T[] result = new T[N];
        foreach (int i in index)
            result[i] = value;
        return result;
    }

    private static int[] findPositiveIndexVector(double[] vector)
    {
        List<int> list = new List<int>();
        for (int i = 0; i < vector.Length; i++)
        {
            if (vector[i] > 0)
            {
                list.Add(i);
            }
        }
        return list.ToArray();
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

    private double tt = 0f;
    private double dt = 0.5f;//0.5f;   // time step
    private double a = 0.12f;      // FHN parameters 
    private double b = 0.013f;     
    private double c1 = 0.26f;      
    private double c2 = 0.1f;
    private double d = 1.0f;
    private double delta = 0.6f;
    private double stimulusStrength = 5f;   // strength of stimulus current
    private double triggerStrength = 3f;   // strength of trigger stimulus current
    private double stimulusTime = 0f;   // time to apply stimulus (in seconds)
    private double[] triggerTime = new double[10];   // time for trigger
    private double stimulusDuration = 10f;   // duration of stimulus (in seconds)
    private int saNodeVertexIndex = 295;   // index of SA node vertex in mesh

    private double[][] laplacian;  // Compute Laplacian operator
    private double[,] lapp;

    private Mesh mesh;   // mesh of heart object
    private Mesh mesh_flip;
    private int N;   // number of vertices in mesh
    private int NF; // number of faces in mesh
    private double[] membranePotentials;   // membrane potentials of vertices
    private double[] gatingVariables;   // gating variables of vertices
    private int[][] triangles;
    private int[][] triangles_flip;
    private Vector3[] vertices;
    private Bounds SANodebounds;
    private Material[] materials_get;
    private Material[] materials_get_flip;
    private MeshFilter meshFilter;
    private MeshFilter meshFilter_flip;
    private MeshRenderer meshRenderer;
    private MeshRenderer meshRenderer_flip;

    private InputData _inputData;
    private bool isParameterEnabled;

    private HashSet<int> ablationNodes;
    private int[] triggerNodes = new int[10];
    private int triggerCurrNum = 0; // trigger number <= 10

    private Matrix<double> m_lap;
    private Vector<double> y;
    private Vector<double> mpv;
    private Vector<double> gvv;
    private Vector<double> I_ext;
    private Vector<double> dydt;
    private Vector<double> y_next;
    private int[] color_idx;

    private GameObject rightHandController;
    private XRRayInteractor xrRayInteractor;
    private Bounds[] vertexbounds;

    private GameObject triggerButtonObject;
    private GameObject triggerTextObject;
    private GameObject leftControllerTextObject;
    private GameObject rightControllerTextObject;
    private Text isAblation;
    private Text leftControllerPosition;
    private Text rightControllerPosition;


    // Start is called before the first frame update
    void Start()
    {
        rightHandController = GameObject.Find("RightHand Controller");
        xrRayInteractor = rightHandController.GetComponent<XRRayInteractor>();
        triggerButtonObject = GameObject.Find("TabletTriggerButton");
        if (triggerButtonObject != null)
        {
            // Get the Button component from the found GameObject
            Button button = triggerButtonObject.GetComponent<Button>();
            button.onClick.AddListener(ButtonClicked);
        }
        triggerTextObject = GameObject.Find("GlassTextForButton");
        if (triggerTextObject != null)
        {
            // Get the Button component from the found GameObject
            isAblation = triggerTextObject.GetComponent<Text>();
        }
        leftControllerTextObject = GameObject.Find("GlassTextForLeftController");
        if (leftControllerTextObject != null)
        {
            // Get the Button component from the found GameObject
            leftControllerPosition = leftControllerTextObject.GetComponent<Text>();
        }
        rightControllerTextObject = GameObject.Find("GlassTextForRightController");
        if (rightControllerTextObject != null)
        {
            // Get the Button component from the found GameObject
            rightControllerPosition = rightControllerTextObject.GetComponent<Text>();
        }
        _inputData = GetComponent<InputData>();
        isParameterEnabled = true;

        // Obj
        GameObject temp = new GameObject();
        temp.name = "parent";
        temp.transform.parent = this.transform;
        materials_get = SetColors(temp);
        triangles = LoadFaces();
        vertices = LoadVertices();

        //SANodebounds = new Bounds(vertices[saNodeVertexIndex], 5);

        mesh = new Mesh();
        mesh.vertices = vertices;
        mesh.subMeshCount = triangles.Length;
        for (int i = 0; i < mesh.subMeshCount; i++)
        {
            mesh.SetTriangles(triangles[i], i);
        }
        N = mesh.vertices.Length;
        NF = mesh.subMeshCount;
        membranePotentials = new double[N];
        gatingVariables = new double[N];

        // MeshFilter
        meshFilter = temp.GetComponent<MeshFilter>();
        if (meshFilter == null) meshFilter = temp.AddComponent<MeshFilter>();
        meshFilter.mesh = mesh;

        // MeshRenderer
        meshRenderer = temp.GetComponent<MeshRenderer>();
        if (meshRenderer == null) meshRenderer = temp.AddComponent<MeshRenderer>();


        // Flip Obj
        GameObject temp_flip = new GameObject();
        temp_flip.name = "flip";
        temp_flip.transform.parent = this.transform;
        materials_get_flip = SetColors(temp_flip);
        triangles_flip = LoadFaces_flip();
        mesh_flip = new Mesh();
        mesh_flip.vertices = vertices;
        mesh_flip.subMeshCount = triangles_flip.Length;
        for (int i = 0; i < mesh_flip.subMeshCount; i++)
        {
            mesh_flip.SetTriangles(triangles_flip[i], i);
        }

        // MeshFilter
        meshFilter_flip = temp_flip.GetComponent<MeshFilter>();
        if (meshFilter_flip == null) meshFilter_flip = temp_flip.AddComponent<MeshFilter>();
        meshFilter_flip.mesh = mesh_flip;

        // MeshRenderer
        meshRenderer_flip = temp_flip.GetComponent<MeshRenderer>();
        if (meshRenderer_flip == null) meshRenderer_flip = temp_flip.AddComponent<MeshRenderer>();

        laplacian = ComputeLaplacian();
        ablationNodes = new HashSet<int>();

        vertexbounds = new Bounds[N];
        for(int i = 0; i < N; i ++)
        {
            vertexbounds[i] = new Bounds(vertices[i], new Vector3(2, 2, 2));
        }
            

        // start
        // Application.targetFrameRate = 60;
    }

    public int[] Advance_realtime(double tt, double dt)
    {
        //laplacian = ComputeLaplacian();

        //Debug.Log("Advance" + tt);

        lapp = new double[N, N];
        for (int i = 0; i < N; i++)
        {
            if (ablationNodes.Contains(i))
                continue;
            int[] k = findNonZeroIndexVector(laplacian[i]);
            foreach(int j in k)
            {
                if (ablationNodes.Contains(j))
                    continue;
                lapp[i, j] = laplacian[i][j];
            }
        }
        m_lap = DenseMatrix.OfArray(lapp);

        // Compute the current state of the membrane potential and gating variables
        y = Vector<double>.Build.Dense(2 * N);
        mpv = Vector<double>.Build.Dense(membranePotentials);
        gvv = Vector<double>.Build.Dense(gatingVariables);
        y.SetSubVector(0, N, mpv);
        y.SetSubVector(N, N, gvv);

        // Compute the external stimulus current
        I_ext = Vector<double>.Build.Dense(N);
        double[] Iex = new double[N];
        int[] idx = new int[N];
        if (tt >= stimulusTime && tt <= stimulusTime + stimulusDuration)
        {
            // Debug.Log("SA node stimulus " + tt + "/" + stimulusDuration);
            idx = findPositiveIndexVector(laplacian[saNodeVertexIndex]); //stimulate SA node every heart cycle 'T'
            //double[] Iex = setSubVector<double>(idx, stimulusStrength);
            //I_ext = Vector<double>.Build.DenseOfArray(Iex);
        }
        foreach (int i in idx)
            Iex[i] = stimulusStrength;
        for(int i = 0; i < triggerCurrNum; i++)
        {
            double trigTime = triggerTime[i];
            if (tt >= trigTime && tt <= trigTime + stimulusDuration)
                idx = findPositiveIndexVector(laplacian[triggerNodes[i]]);
        }
        foreach (int i in idx)
            Iex[i] = triggerStrength;
        I_ext = Vector<double>.Build.DenseOfArray(Iex);
        dydt = ComputeDerivatives(y, m_lap, I_ext);

        // Integrate the state variables using Forward Euler method
        y_next = y + dt * dydt;

        // Update the membrane potentials and gating variables
        membranePotentials = y_next.SubVector(0, N).ToArray();
        gatingVariables = y_next.SubVector(N, N).ToArray();

        // Display whole heart
        Vector<double> u = Vector<double>.Build.DenseOfArray(membranePotentials);
        Vector<double> m = 1 + roundVector(63 * u);
        Vector<double> ones = Vector<double>.Build.Dense(N, 1);
        Vector<double> sixtyfour = Vector<double>.Build.Dense(N, 64);
        m = m.PointwiseMaximum(ones); // m = max(m, 1); 
        m = m.PointwiseMinimum(sixtyfour); // min(m, 64);
        double[] v_colors = (m - 1).ToArray();
        int[] color_idx = v_colors.Select(d => (int)d).ToArray(); // convert a double array to integer one
        return color_idx;
    }

    private void ButtonClicked()
    {
        Debug.Log("Button pressed! - " + isParameterEnabled);
        isParameterEnabled ^= true;
    }

    void Update()
    {
        // laplacian = ComputeLaplacian();
        if (tt < 4800)
        {
            int[] color_idx = Advance_realtime(tt, dt);
            
            // Update the current time
            tt += dt;

            // ablation nodes visualization
            foreach (int i in ablationNodes)
                color_idx[i] = 55;

            // trigger nodes visualization
            foreach (int i in triggerNodes)
                color_idx[i] = 18;

            // visualize
            ShowGameObject(triangles, mesh, meshRenderer, materials_get, color_idx);
            ShowGameObject(triangles, mesh_flip, meshRenderer_flip, materials_get_flip, color_idx);
        }

        /*if (Input.GetKeyDown(KeyCode.Space))
        {
            isParameterEnabled ^= true;
        }*/

        if (isParameterEnabled)
        {
            isAblation.text = "Ablation Mode";
            Ray ray = new Ray(xrRayInteractor.transform.position, xrRayInteractor.transform.forward);
            for (int i = 0; i < N; i++)
            {
                if(vertexbounds[i].IntersectRay(ray, out float intersectionDistance))
                {
                    if (intersectionDistance < 10f && ablationNodes.Count <= 100)
                    {
                        Debug.Log("intersection distance " + intersectionDistance);
                        ablationNodes.Add(i);
                        Debug.Log(ablationNodes.Count + "/100, touching ablation node " + i);
                    }
                }
            }
        }
        else
        {
            isAblation.text = "Trigger Mode";
            Ray ray = new Ray(xrRayInteractor.transform.position, xrRayInteractor.transform.forward);
            for (int i = 0; i < N; i++)
            {
                if (vertexbounds[i].IntersectRay(ray, out float intersectionDistance))
                {
                    if (intersectionDistance < 10f && triggerCurrNum < 10)
                    {
                        Debug.Log("intersection distance " + intersectionDistance);
                        triggerNodes[triggerCurrNum] = i;
                        triggerTime[triggerCurrNum] = tt;
                        triggerCurrNum += 1;
                        Debug.Log(triggerCurrNum + "/10, node " + i + " as trigger");
                    }
                }
            }
        }


        if (_inputData._rightController.TryGetFeatureValue(CommonUsages.devicePosition, out Vector3 rightPosition))
        {
            rightControllerPosition.text = "Right Controller Position: "+ rightPosition.ToString("F2");
        }
        if (_inputData._leftController.TryGetFeatureValue(CommonUsages.devicePosition, out Vector3 leftPosition))
        {
            leftControllerPosition.text = "Left Controller Position: " + leftPosition.ToString("F2");
        }
    }


    private static Vector<double> roundVector(Vector<double> vector)
    {
        double[] arr_src = vector.ToArray();
        double[] arr_dst = new double[vector.Count];
        for (int i = 0; i < vector.Count; i++)
        {
            arr_dst[i] = Math.Round(arr_src[i], MidpointRounding.AwayFromZero);
        }
        Vector<double> result = Vector<double>.Build.DenseOfArray(arr_dst); ;
        return result;
    }

    private void ShowGameObject(int[][] faces, Mesh mesh, MeshRenderer meshRenderer, Material[] materials_get, int[] verticesColor)
    {
        // Debug.Log("ShowGameObject");
        int[] color_idx = ColorV2F(verticesColor, faces);
        Material[] materials_push = new Material[mesh.subMeshCount];
        for (int i = 0; i < mesh.subMeshCount; i++)
        {
            materials_push[i] = materials_get[color_idx[i]];
        }
        meshRenderer.materials = materials_push;
    }


    //Matrix<double> ComputeLaplacian()
    double[][] ComputeLaplacian()
    {
        double[][] edge = new double[N][];
        for (int i = 0; i < N; i++)
        {
            edge[i] = new double[N];
        }
        for (int i = 0; i < NF; i++)
        {
            Matrix A = DenseMatrix.OfArray(new double[,] {
                { vertices[triangles[i][0]][0], vertices[triangles[i][0]][1], vertices[triangles[i][0]][2]},
                { vertices[triangles[i][1]][0], vertices[triangles[i][1]][1], vertices[triangles[i][1]][2]},
                { vertices[triangles[i][2]][0], vertices[triangles[i][2]][1], vertices[triangles[i][2]][2]}});

            Matrix B = DenseMatrix.OfArray(new double[,] {
                { vertices[triangles[i][1]][0], vertices[triangles[i][1]][1], vertices[triangles[i][1]][2]},
                { vertices[triangles[i][2]][0], vertices[triangles[i][2]][1], vertices[triangles[i][2]][2]},
                { vertices[triangles[i][0]][0], vertices[triangles[i][0]][1], vertices[triangles[i][0]][2]}});

            var Diff = A.Subtract(B);
            var Norm = Diff.PointwisePower(2).RowSums().PointwiseSqrt();

            edge[triangles[i][0]][triangles[i][1]] = Norm[0];
            edge[triangles[i][1]][triangles[i][2]] = Norm[1];
            edge[triangles[i][2]][triangles[i][0]] = Norm[2];

            edge[triangles[i][1]][triangles[i][0]] = Norm[0];
            edge[triangles[i][2]][triangles[i][1]] = Norm[1];
            edge[triangles[i][0]][triangles[i][2]] = Norm[2];
        }

        double[][] lap = new double[N][];
        for (int i = 0; i < N; i++)
        {
            lap[i] = new double[N];
        }
        for (int i = 0; i < N; i++)
        {
            int[] k = findNonZeroIndexVector(edge[i]);
            int ni = k.Length; // the number of neighbours
            if (ni != 0)
            {
                double hi = avgNeighbourDist(edge[i], k, ni); // the average distance to the neighbours
                double invhi = avgNeighbourDistInverse(edge[i], k, ni); // the average inverse distance to the neighbours

                // Laplacian of vertex itself
                lap[i][i] = -(4 / hi) * invhi;

                // Laplacian of direct neighbours
                foreach (int j in k)
                {
                    lap[i][j] = 1.0 / edge[i][j] * (4 / (hi * ni));
                }
            }
        }
        return lap;
    }

    Vector<double> ComputeDerivatives(Vector<double> y, Matrix<double> laplacian, Vector<double> I_ext)
    {
        Vector<double> dydt = Vector<double>.Build.Dense(2 * N);

        Vector<double> V = y.SubVector(0, N);
        Vector<double> W = y.SubVector(N, N);

        Vector<double> dVdt = Vector<double>.Build.Dense(N);
        Vector<double> dWdt = Vector<double>.Build.Dense(N);

        // Compute dV/dt
        dVdt = c1 * V.PointwiseMultiply(V - a).PointwiseMultiply(1 - V) - c2 * V.PointwiseMultiply(W) + I_ext + delta * laplacian * V;
        // Debug.Log("dVm" + dVdt[saNodeVertexIndex]);

        // Compute dW/dt
        dWdt = b * (V - d * W);
        // Ablation
        /*
        for(int i = 0; i < 100; i++)
        {
            dWdt.SetSubVector(saNodeVertexIndex + i, 1, Vector<double>.Build.DenseOfArray(new double[] { 0f }));
        }
        Debug.Log("dW" + dWdt[saNodeVertexIndex + 50]);*/

        dydt.SetSubVector(0, N, dVdt);
        dydt.SetSubVector(N, N, dWdt);

        return dydt;
    }
}


