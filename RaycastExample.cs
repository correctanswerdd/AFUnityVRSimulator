using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.XR;
using UnityEngine.XR.Interaction.Toolkit;

public class RaycastExample : MonoBehaviour
{
    private XRRayInteractor xrrayinteractor;
    private Vector3[] vertices;
    private int saNodeVertexIndex = 295;   // index of SA node vertex in mesh
    private Bounds SANodebounds;

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

    // Start is called before the first frame update
    void Start()
    {
        xrrayinteractor = GetComponent<XRRayInteractor>();
        vertices = LoadVertices();
        SANodebounds = new Bounds(vertices[saNodeVertexIndex], new Vector3(5, 5, 5));
    }

    // Update is called once per frame
    void Update()
    {
        /*RaycastHit res;
        if (xrrayinteractor.TryGetCurrent3DRaycastHit(out res))
        {
            Vector3 groundPt = res.point; // the coordinate that the ray hits
            //Debug.Log(" coordinates on the ground: " + groundPt);
        }*/

        Ray ray = new Ray(xrrayinteractor.transform.position, xrrayinteractor.transform.forward);
        if (SANodebounds.IntersectRay(ray))
        {
            Debug.Log("touch SA node");
        }

    }
}
