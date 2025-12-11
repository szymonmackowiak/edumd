

BufferedReader reader1, reader2, reader3;
String lineX, lineY, lineZ;
 
int xsize=1000, ysize=1000, choice;
float skala=200;
float x0, y0, z0, eyeX, eyeY, eyeZ;
float teta, fi, r, deltaAlfa, deltaR, sphR;
float teta0, fi0, r0, sphR0;

void setup() {
  // Open the file from the createWriter() example
  choice=1;
  reader1 = createReader("../output/RX.txt");  
  reader2 = createReader("../output/RY.txt"); 
  reader3 = createReader("../output/RZ.txt"); 
  frameRate(20);
  size(600,600,P3D);
  background(255);  
  //camera(100, 100, 200, 100, 100, 200, 1, 1, 1);
  r0=500;
  teta0=0;
  fi0=0;
  sphR0=5;
  
  if(choice==2)
  {
  teta0=6.28*0.25;
  r0=50;
  }
  
  r=r0;
  teta=teta0;
  fi=fi0;
  sphR=sphR0;
  
  x0=0.5*xsize;
  y0=0.5*ysize;
  z0=0.5*xsize;
  deltaAlfa=0.05;
  deltaR=50;
  
}
void draw() 
{
  
  try 
  {
    lineX = reader1.readLine();
    lineY = reader2.readLine();
    lineZ = reader3.readLine();
  } 
  catch (IOException e) 
  {
    e.printStackTrace();
    lineX = null;
    lineY = null;
    lineZ = null;
  }
  if (lineX == null || lineY == null || lineZ == null) 
  {
    //noLoop();
  setup();  
  } 
  else 
  {
    String[] piecesX = split(lineX, " ");
    String[] piecesY = split(lineY, " ");
    String[] piecesZ = split(lineZ, " ");
    
    float[] x= new float[108];
    float[] y= new float[108];
    float[] z= new float[108];
    
    
    for(int i=0; i<108; i++)
    {
    x[i] = float(piecesX[i])*20+500;
    y[i] = float(piecesY[i])*20+500;
    z[i] = float(piecesZ[i])*20+500;
    }
    
    
    background(255);  
    noStroke();
    lights();

if (keyPressed) 
  {
    if (key == 'a') 
    {
      teta=teta-deltaAlfa;
    }
    else if (key == 'd') 
    {
      teta=teta+deltaAlfa;
    }
    else if (key == 'w') 
    {
      fi=fi-deltaAlfa;
    }
    else if (key == 's') 
    {
      fi=fi+deltaAlfa;
    }
    else if (key == 'q') 
    {
      r=r+deltaR;
    }
    else if (key == 'e') 
    {
      r=r-deltaR;
    }  
    else if (key == 'r') 
    {
      sphR=sphR+0.5;
    }
    else if (key == 'f') 
    {
      sphR=sphR-0.5;
    }
    else if (key == 'c') 
    {
      r=r0;
      teta=teta0;
      fi=fi0;
      sphR=sphR0;
    }
    
  } 
  else 
  {
    teta=teta;
    fi=fi;
    r=r;
  }



eyeX=x0+r*cos(teta)*cos(fi);
eyeY=y0+r*cos(teta)*sin(fi);
eyeZ=z0+r*sin(teta);


if(choice==1)
{
camera(eyeX, eyeY, eyeZ, // eyeX, eyeY, eyeZ
         0.5*xsize, 0.5*ysize, 0.5*ysize, // centerX, centerY, centerZ
         0.0, 1.0, 0.0); // upX, upY, upZ
}
    //ellipse(x1,y1,skala,skala);
    //ellipse(x2,y2,skala,skala);
    //ellipse(x3,y3,400,400);
    
if(choice==3)
{
ortho(-100, 500, -500, 100);    
}
    
if(choice==1)
{
    
   
    
    
    for(int i=0; i<108; i++)
    {
      if (i==0){
      fill(250, 20, 50);
      }
      else{
        fill(50, 205, 50);
      }
      
    pushMatrix();
    translate(x[i], z[i], y[i]);
    sphereDetail(10);
    sphere(sphR);
    popMatrix();
    }
    
    
}
    
else if(choice==2)
{
    teta=6.28/4;
    
    
    fill(50, 205, 50);
    for(int i=0; i<108; i++)
    {
    ellipse(y[i],z[i],sphR,sphR);
    }
    
    
}

else if(choice==3)
{
    
    
    fill(50, 205, 50);
    for(int i=0; i<108; i++)
    {
    pushMatrix();
    translate(y[i], z[i], x[i]);
    sphereDetail(10);
    sphere(sphR);
    popMatrix();
    }
    
    
}
    
    
  }
} 
