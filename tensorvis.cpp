/**********************************************************************
  TensorVis -- generate and display tensor representations

  Copyright (C) 2011 David C. Lonie

  This source code is released under the New BSD License, (the "License").

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

 **********************************************************************/

#include "tensorvis.h"

#include "ui/tensorviswidget.h"

#include <avogadro/color3f.h>
#include <avogadro/cube.h>
#include <avogadro/glwidget.h>
#include <avogadro/mesh.h>
#include <avogadro/molecule.h>
#include <avogadro/painter.h>

#include <QtCore/QDebug>
#include <QtCore/QTime>
#include <QtGui/QAction>
#include <QtGui/QDockWidget>

#include <Eigen/Geometry>

#include <float.h>

#undef DEBUG_TENSORGEN

using namespace Avogadro;

namespace TensorVis {

  TensorVis::TensorVis(QObject *parent)
    : m_molecule(0),
      m_dock(0)
  {
    QAction *action = new QAction(this);
    action->setText(tr("Add &Tensor..."));
    m_actions.append(action);
  }

  TensorVis::~TensorVis()
  {
  }

  QList<QAction *> TensorVis::actions() const
  {
    return m_actions;
  }

  QString TensorVis::menuPath(QAction *action) const
  {
    return tr("E&xtensions");
  }

  // "point" is on a sphere formed by iterating over phi and theta
  // in spherical coordinates:
  // x = r * sin(theta)*cos(phi)
  // y = r * sin(theta)*sin(phi)
  // z = r * cos(theta)
  inline float TensorVis::calculateTensorValue(const Eigen::Vector3f &point,
                                               const Eigen::Matrix3f &tensor)
  {
    return (point.transpose() * tensor).dot(point);
  }

  QUndoCommand * TensorVis::performAction(QAction *action,
                                          GLWidget *widget)
  {
    if (!m_dock) {
      qDebug() << "Dock widget not yet created?";
      return 0;
    }
    m_dock->show();
    return 0;
  }

  void TensorVis::generateMesh()
  {
    if (!m_dock) {
      qDebug() << "Dock widget not yet created?";
      return;
    }

    QTime startTime (QTime::currentTime());

    const Eigen::Matrix3f tensor (m_dock->tensor());
    const Eigen::Vector3f tensorOrigin (m_dock->origin());

    const Color3f posColor (m_dock->posColor());
    const Color3f negColor (m_dock->negColor());

    float THETA_STEP, PHI_STEP;
    switch(this->m_dock->resolution())
    {
    case Poor:
      THETA_STEP = PHI_STEP = 0.25;
      break;
    default:
    case Low:
      THETA_STEP = PHI_STEP = 0.15;
      break;
    case Medium:
      THETA_STEP = PHI_STEP = 0.10;
      break;
    case High:
      THETA_STEP = PHI_STEP = 0.05;
      break;
    case Ultra:
      THETA_STEP = PHI_STEP = 0.01;
      break;
    }

    const float TENSOR_RADIUS = this->m_dock->scale();

    // Theta runs [0,pi] (note that pi is included)
    const unsigned long NUM_THETA =
      static_cast<unsigned long>(M_PI / THETA_STEP) + 1;
    // Phi runs [0,2*pi) (note that 2*pi is not included)
    const unsigned long NUM_PHI =
      static_cast<unsigned long>(2*M_PI / PHI_STEP);
    const unsigned long NUM_POINTS = NUM_THETA * NUM_PHI;

    // Set up caches
    QVector<float> sinPhis;
    sinPhis.reserve(NUM_PHI);
    QVector<float> cosPhis;
    cosPhis.reserve(NUM_PHI);
    QVector<float> sinThetas;
    sinThetas.reserve(NUM_THETA);
    QVector<float> cosThetas;
    cosThetas.reserve(NUM_THETA);
    QVector<Eigen::Vector3f> tensorPoints (NUM_POINTS);
    QVector<float> tensorValue (NUM_POINTS);

    // Populate caches
    // Subtract 1/2 step to drown out floating point noise
    const float thetaMax = ((NUM_THETA - 0.5) * THETA_STEP);
    for (float theta = 0.0; theta < thetaMax; theta += THETA_STEP) {
      sinThetas.push_back(sin(theta));
      cosThetas.push_back(cos(theta));
    }
    // Subtract 1/2 step to drown out floating point noise
    const float phiMax = ((NUM_PHI - 0.5) * PHI_STEP);
    for (float phi = 0.0; phi < phiMax; phi += PHI_STEP) {
      sinPhis.push_back(sin(phi));
      cosPhis.push_back(cos(phi));
    }

    // ensure that the first and last entries of the lists are the
    // same (sin and cos can behave strangely on different
    // implementations). Only do this for the theta values, as they
    // include both endpoints.
    sinThetas.last() = sinThetas.first();
    cosThetas.last() = -cosThetas.first();

    // Check that all spherical coordinates were cached
#ifdef DEBUG_TENSORGEN
    printf("sinThetas: %u num: %u\n", sinThetas.size(), NUM_THETA);
    printf("sinPhis:   %u num: %u\n", sinPhis.size(), NUM_PHI);
#endif
    Q_ASSERT(sinPhis.size() == NUM_PHI);
    Q_ASSERT(sinThetas.size() == NUM_THETA);

    unsigned long pointInd = 0;
    Eigen::Vector3f tmpVec;
    double value;
    for (unsigned long thetaInd = 0; thetaInd < NUM_THETA; ++thetaInd) {
      const float &cosTheta = cosThetas[thetaInd];
      const float &sinTheta = sinThetas[thetaInd];
      for (unsigned long phiInd = 0; phiInd < NUM_PHI; ++phiInd) {
        const float &cosPhi = cosPhis[phiInd];
        const float &sinPhi = sinPhis[phiInd];
        tmpVec = Eigen::Vector3f(sinTheta*cosPhi,
                                 sinTheta*sinPhi,
                                 cosTheta);
        value = calculateTensorValue(tmpVec, tensor);
        tmpVec *= value;
        tmpVec *= TENSOR_RADIUS;

        tensorPoints[pointInd] = tmpVec + tensorOrigin;
        tensorValue[pointInd] = value;
        ++pointInd;
      }
    }

    // Check that all points were generated
    Q_ASSERT(tensorPoints.size() == NUM_POINTS);

    // Generate mesh -- expects triplets of vertices of clockwise
    // wound triangles
    //
    //     ^
    //  t  |
    //  h  |  2--4
    //  e  |  |\ |
    //  t  |  | \|
    //  a  |  1--3
    //     |
    //     +------->
    //        phi
    //
    //  Iterating over point with counter ind,
    //
    //  p1 = points[ind]
    //  p2 = points[ind + NUM_PHI]
    //  p3 = points[ind + 1]
    //  p4 = points[ind + 1 + NUM_PHI]
    //
    // Triangles (p1, p2, p3) and (p2, p4, p3) will cover this space.

    std::vector<Eigen::Vector3f> verts;
    std::vector<Eigen::Vector3f> norms;
    std::vector<Color3f> colors;
    // @todo reserve appropriate amount of space

    // This vector is used to smooth out the normals of the surface. It
    // stores an index to each normal added to norms. The outer vector is
    // indexed by tensorPoints index, and the inner vectors contain all
    // normals calculated at that point. After calculating the vertices and
    // normals, the Eigen::Vectors in each inner vector of normalSmoother are
    // averaged and set to the average. The improves the appearance of the
    // final mesh.
    std::vector<std::vector<size_t> > normalSmoother (NUM_POINTS);

    // Points:
    const Eigen::Vector3f *p1, *p2, *p3, *p4;
    // Normals:
    Eigen::Vector3f quadNormal;
    // Temp vecs used in normal calculations:
    Eigen::Vector3f tv13, tv14, tv23, tv24;
    // Values:
    const float *v1, *v2, *v3, *v4;
    // Indices:
    unsigned long i1, i2, i3, i4;
    bool onlyRenderOneTriangle;

    // When indicies are above this threshold they must be rendered
    // differently:
    const unsigned long autoGenMeshLimit = NUM_POINTS - NUM_PHI - 1;

    // Build the mesh:
    const unsigned long meshCalcStart = 0;
    const unsigned long meshCalcEnd = NUM_POINTS - NUM_PHI;
    for (unsigned long ind = meshCalcStart; ind < meshCalcEnd; ++ind) {
      // Initialize
      onlyRenderOneTriangle = false;

      // Assign indices
      i1 = ind;
      i2 = ind + NUM_PHI;
      i3 = ind + 1;
      i4 = ind + NUM_PHI + 1;

      // If we are on a boundary, subtract NUM_PHI from indices 3 and 4.
      if ((ind + 1) % NUM_PHI == 0) {
        i3 -= NUM_PHI;
        i4 -= NUM_PHI;
      }

      // Adjust the last loop so that we don't exceed tensorPoint's
      // bounds. Just set i2 and i4 to NUM_POINTS-1, since the last
      // row of points all have the same coords.
      if (ind >= autoGenMeshLimit) {
        // For this last loop, we'll only render a single triangle at a time
        onlyRenderOneTriangle = true;
        // Skip the last few points since there's not enough data to
        // draw the quad.
        if (ind >= NUM_POINTS - 4) break;
        // i4 will not be used -- set to a safe index
        i4 = 0;
      }

      // Assign points
      p1 = &tensorPoints[i1];
      p2 = &tensorPoints[i2];
      p3 = &tensorPoints[i3];
      p4 = &tensorPoints[i4];

#ifdef DEBUG_TENSORGEN
      printf("\nNew quad:\n");
      printf("  %9.5f %9.5f %9.5f\n", p1->x(), p1->y(), p1->z());
      printf("  %9.5f %9.5f %9.5f\n", p2->x(), p2->y(), p2->z());
      printf("  %9.5f %9.5f %9.5f\n", p3->x(), p3->y(), p3->z());
      printf("  %9.5f %9.5f %9.5f\n", p4->x(), p4->y(), p4->z());
#endif

      // Assign normals.
      tv13 = *p1 - *p3;
      tv14 = *p1 - *p4;
      tv23 = *p2 - *p3;
      tv24 = *p2 - *p4;
      bool isZero13 = tv13.isZero(1e-5);
      bool isZero14 = tv14.isZero(1e-5);
      bool isZero23 = tv23.isZero(1e-5);
      bool isZero24 = tv24.isZero(1e-5);
      // Normalizations are performed in the later smoothing step. Note that
      // we only include p4 if !onlyRenderOneTriangle
      quadNormal << 0.0, 0.0, 0.0;
      if (!onlyRenderOneTriangle) {
        if (!isZero14 && !isZero13)
          quadNormal += tv14.cross(tv13);
        if (!isZero24 && !isZero23)
          quadNormal += tv24.cross(tv23);
        if (!isZero13 && !isZero23)
          quadNormal += tv13.cross(tv23);
        if (!isZero14 && !isZero24)
          quadNormal += tv14.cross(tv24);
      }
      else {
        if (!isZero13 && !isZero23)
          quadNormal += tv13.cross(tv23);
      }
      if (quadNormal.isZero(1e-5))
        quadNormal = *p1 - tensorOrigin; // All points are the same. Approx.

#ifdef DEBUG_TENSORGEN
      printf("\nQuadNormal: %9.5f %9.5f %9.5f\n",
             quadNormal.x(), quadNormal.y(), quadNormal.z());
#endif

      // Assign values
      v1 = &tensorValue[i1];
      v2 = &tensorValue[i2];
      v3 = &tensorValue[i3];
      v4 = &tensorValue[i4];

      // If the shared vertex is negative, reverse the winding:
      if (*v2 >= 0 && *v3 >= 0) {
        verts.push_back(*p1);
        verts.push_back(*p2);
        verts.push_back(*p3);
        norms.push_back(quadNormal);
        normalSmoother[i1].push_back(norms.size() - 1);
        norms.push_back(quadNormal);
        normalSmoother[i2].push_back(norms.size() - 1);
        norms.push_back(quadNormal);
        normalSmoother[i3].push_back(norms.size() - 1);
        colors.push_back(posColor);
        colors.push_back(posColor);
        colors.push_back(posColor);

        if (!onlyRenderOneTriangle) {
          verts.push_back(*p2);
          verts.push_back(*p4);
          verts.push_back(*p3);
          norms.push_back(quadNormal);
          normalSmoother[i2].push_back(norms.size() - 1);
          norms.push_back(quadNormal);
          normalSmoother[i4].push_back(norms.size() - 1);
          norms.push_back(quadNormal);
          normalSmoother[i3].push_back(norms.size() - 1);
          colors.push_back(posColor);
          colors.push_back(posColor);
          colors.push_back(posColor);
        }
      }
      else {
        verts.push_back(*p3);
        verts.push_back(*p2);
        verts.push_back(*p1);
        norms.push_back(-quadNormal);
        normalSmoother[i3].push_back(norms.size() - 1);
        norms.push_back(-quadNormal);
        normalSmoother[i2].push_back(norms.size() - 1);
        norms.push_back(-quadNormal);
        normalSmoother[i1].push_back(norms.size() - 1);
        colors.push_back(negColor);
        colors.push_back(negColor);
        colors.push_back(negColor);

        if (!onlyRenderOneTriangle) {
          verts.push_back(*p3);
          verts.push_back(*p4);
          verts.push_back(*p2);
          norms.push_back(-quadNormal);
          normalSmoother[i3].push_back(norms.size() - 1);
          norms.push_back(-quadNormal);
          normalSmoother[i4].push_back(norms.size() - 1);
          norms.push_back(-quadNormal);
          normalSmoother[i2].push_back(norms.size() - 1);
          colors.push_back(negColor);
          colors.push_back(negColor);
          colors.push_back(negColor);
        }
      }
    }

    // Smooth out the normals:
    for (std::vector<std::vector<size_t> >::const_iterator
         it_out = normalSmoother.begin(), it_out_end = normalSmoother.end();
         it_out != it_out_end; ++it_out) {
      tmpVec << 0.0,0.0,0.0;
#ifdef DEBUG_TENSORGEN
      printf("\nNew set of vectors to smooth!\n");
#endif
      // Calculate average. Note that we don't keep track of how many vectors
      // are being averaged (i.e. the denominator in the mean calculation).
      // Since we are normalizing the averaged vector anyway, the division
      // would be a waste of cycles.
      for (std::vector<size_t>::const_iterator it = it_out->begin(),
           it_end = it_out->end(); it != it_end; ++it) {
        tmpVec += norms[*it];
#ifdef DEBUG_TENSORGEN
        printf("  Vector #%5d: %9.5f %9.5f %9.5f\n", *it,
               norms[*it].x(), norms[*it].y(), norms[*it].z());
#endif
      }
      tmpVec.normalize();
#ifdef DEBUG_TENSORGEN
      printf("  Average:   %9.5f %9.5f %9.5f Length: %9.5f\n",
             tmpVec.x(), tmpVec.y(), tmpVec.z(), tmpVec.norm());
#endif
      // Reset the vectors
      for (std::vector<size_t>::const_iterator it = it_out->begin(),
           it_end = it_out->end(); it != it_end; ++it) {
        norms[*it] = tmpVec;
      }
    }

    qDebug() << "Tensor mesh generated:" << norms.size() / 3 << "triangles,"
             << verts.size() << "vertices," << norms.size() << "normals."
             << "Time elapsed (s):"
             << startTime.msecsTo(QTime::currentTime()) / 1000.0;

    // Look for existing Tensor mesh in molecule
    QList<Mesh*> meshes = this->m_molecule->meshes();
    Mesh *mesh = NULL;
    for (QList<Mesh*>::const_iterator it = meshes.constBegin(),
         it_end = meshes.constEnd(); it != it_end; ++it) {
      if ((*it)->name().startsWith("Tensor")) {
        mesh = *it;
        break;
      }
    }
    if (mesh == NULL)
      mesh = m_molecule->addMesh();

    mesh->setVertices(verts);
    mesh->setNormals(norms);
    mesh->setColors(colors);
    mesh->setName(tr("Tensor at (%L1, %L2, %L3), trace average = %L4")
                  .arg(tensorOrigin.x(), 1, 'f')
                  .arg(tensorOrigin.y(), 1, 'f')
                  .arg(tensorOrigin.z(), 1, 'f')
                  .arg(tensor.trace() * 0.33333333));

    // Work around a bug in surfaceengine: All meshes must have cubes
    // or that engine can crash.
    if (!mesh->cube()) {
      Cube *cube = m_molecule->addCube();
      cube->setLimits(Eigen::Vector3d(0,0,0), Eigen::Vector3i(1,1,1), 1);
      std::vector<double> junkDVec (1);
      cube->setData(junkDVec);
      mesh->setCube(cube->id());
    }

    this->m_molecule->update();
  }

  void TensorVis::clearMesh()
  {
    // Look for existing Tensor mesh in molecule
    QList<Mesh*> meshes = this->m_molecule->meshes();
    Mesh *mesh = NULL;
    for (QList<Mesh*>::const_iterator it = meshes.constBegin(),
         it_end = meshes.constEnd(); it != it_end; ++it) {
      if ((*it)->name().startsWith("Tensor")) {
        mesh = *it;
        break;
      }
    }
    if (mesh != NULL)
      this->m_molecule->removeMesh(mesh);
    this->m_molecule->update();
  }

  void TensorVis::setMolecule(Molecule *molecule)
  {
    m_molecule = molecule;
  }

  QDockWidget * TensorVis::dockWidget()
  {
    if (!m_dock) {
      m_dock = new TensorVisWidget ();
      connect(m_dock, SIGNAL(generateMesh()), SLOT(generateMesh()));
      connect(m_dock, SIGNAL(clearMesh()), SLOT(clearMesh()));
    }
    m_dock->hide();
    return m_dock;
  }

}

#include "tensorvis.moc"

Q_EXPORT_PLUGIN2(tensorvis, TensorVis::TensorVisFactory)
