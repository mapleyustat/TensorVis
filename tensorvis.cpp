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
  float TensorVis::calculateTensorValue(const Eigen::Vector3f &point,
                                        const Eigen::Matrix3f &tensor)
  {
    return (point.transpose() * tensor).dot(point);
  }

  QUndoCommand * TensorVis::performAction(QAction *action,
                                          GLWidget *widget)
  {
    QTime startTime (QTime::currentTime());
    /// @todo these should be set by user eventually:
    Eigen::Matrix3f tensor;
    tensor <<
      88.921500, 46.045900, -0.777100,
      17.744000, 18.222400, -0.248500,
      -0.900100, -0.432300, -81.653700;
    const Eigen::Vector3f tensorOrigin (0.0, 0.0, 0.0);

    const float posR = 0.2;
    const float posG = 0.2;
    const float posB = 0.6;
    const float negR = 0.6;
    const float negG = 0.2;
    const float negB = 0.2;
    /// @todo Predefine resolutions
    const float THETA_STEP = 0.05;
    const float PHI_STEP = 0.05;
    const float TENSOR_RADIUS = .05;

    const float TWO_PI = 2*M_PI;
    const unsigned long NUM_THETA =
      static_cast<unsigned long>(M_PI / THETA_STEP);
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
    QVector<Eigen::Vector3f> tensorNormals (NUM_POINTS);
    QVector<float> tensorValue (NUM_POINTS);

    // Populate caches
    /// @todo exploit symmetry?
    for (float theta = 0.0; theta <= M_PI; theta += THETA_STEP) {
      sinThetas.push_back(sin(theta));
      cosThetas.push_back(cos(theta));
    }
    for (float phi = 0.0; phi <= TWO_PI; phi += PHI_STEP) {
      sinPhis.push_back(sin(phi));
      cosPhis.push_back(cos(phi));
    }

    // ensure that the first and last entries of the phi lists are the
    // same (sin and cos can behave strangely on different
    // implementations)
    sinPhis.first() = sinPhis.last();
    cosPhis.first() = cosPhis.last();

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

        // pointInd will overflow to 0 on first iteration (prefix
        // increment is faster)
        tensorPoints[pointInd] = tmpVec + tensorOrigin;
        tensorNormals[pointInd] = tmpVec - tensorOrigin;
        tensorValue[pointInd] = value;
        ++pointInd;
      }
    }

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

    const Eigen::Vector3f *p1, *p2, *p3, *p4;
    const Eigen::Vector3f *n1, *n2, *n3, *n4;
    const float *v1, *v2, *v3, *v4;
    unsigned long i1, i2, i3, i4;

    // When indicies are above this threshold they must be rendered
    // differently:
    const unsigned long autoGenMeshLimit = NUM_POINTS - NUM_PHI - 1;

    // Build the mesh:
    for (unsigned long ind = 0; ind < NUM_POINTS; ++ind) {
      // Assign indices
      i1 = ind;
      i2 = ind + NUM_PHI;
      i3 = ind + 1;
      i4 = ind + NUM_PHI + 1;

      // If we are on a bondary, subtract NUM_PHI from indices 3 and 4.
      if (ind != 0 && (ind % NUM_PHI == 0)) {
        i3 -= NUM_PHI;
        i4 -= NUM_PHI;
      }

      // Adjust the last loop so that we don't exceed tensorPoint's
      // bounds. Just set i2 and i4 to NUM_POINTS-1, since the last
      // row of points all have the same coords.
      if (ind >= autoGenMeshLimit) {
        // Skip the last few points since there's not enough data to
        // draw the quad.
        if (ind >= NUM_POINTS - 4) break;
        i2 = NUM_POINTS-1;
        i4 = NUM_POINTS-1;
      }

      // Assign points
      p1 = &tensorPoints[i1];
      p2 = &tensorPoints[i2];
      p3 = &tensorPoints[i3];
      p4 = &tensorPoints[i4];

      // Assign normals
      n1 = &tensorNormals[i1];
      n2 = &tensorNormals[i2];
      n3 = &tensorNormals[i3];
      n4 = &tensorNormals[i4];

      // Assign values
      v1 = &tensorValue[i1];
      v2 = &tensorValue[i2];
      v3 = &tensorValue[i3];
      v4 = &tensorValue[i4];

      // If the normal shared vertex is negative, reverse the winding:
      if (*v2 >= 0 && *v3 >= 0) {
        verts.push_back(*p1);
        verts.push_back(*p2);
        verts.push_back(*p3);
        norms.push_back(*n1);
        norms.push_back(*n2);
        norms.push_back(*n3);
        /// @todo Scale colors?
        colors.push_back(Color3f(posR, posG, posB));
        colors.push_back(Color3f(posR, posG, posB));
        colors.push_back(Color3f(posR, posG, posB));

        verts.push_back(*p2);
        verts.push_back(*p4);
        verts.push_back(*p3);
        norms.push_back(*n2);
        norms.push_back(*n4);
        norms.push_back(*n3);
        colors.push_back(Color3f(posR, posG, posB));
        colors.push_back(Color3f(posR, posG, posB));
        colors.push_back(Color3f(posR, posG, posB));
      }
      else {
        verts.push_back(*p3);
        verts.push_back(*p2);
        verts.push_back(*p1);
        norms.push_back(*n3);
        norms.push_back(*n2);
        norms.push_back(*n1);
        colors.push_back(Color3f(negR, negG, negB));
        colors.push_back(Color3f(negR, negG, negB));
        colors.push_back(Color3f(negR, negG, negB));

        verts.push_back(*p3);
        verts.push_back(*p4);
        verts.push_back(*p2);
        norms.push_back(*n3);
        norms.push_back(*n4);
        norms.push_back(*n2);
        colors.push_back(Color3f(negR, negG, negB));
        colors.push_back(Color3f(negR, negG, negB));
        colors.push_back(Color3f(negR, negG, negB));
      }
    }

    qDebug() << "Tensor mesh generated:" << norms.size() / 3 << "triangles,"
             << verts.size() << "vertices," << norms.size() << "normals."
             << "Time elapsed (s):"
             << startTime.msecsTo(QTime::currentTime()) / 1000.0;

    Mesh *mesh = m_molecule->addMesh();
    mesh->setVertices(verts);
    mesh->setNormals(norms);
    mesh->setColors(colors);
    /// @todo Ask Jochen if there is an easier way for users to
    /// identify the tensors.
    mesh->setName(tr("Tensor at (%L1, %L2, %L3), trace average = %L4")
                  .arg(tensorOrigin.x(), 1, 'f')
                  .arg(tensorOrigin.y(), 1, 'f')
                  .arg(tensorOrigin.z(), 1, 'f')
                  .arg(tensor.trace() * 0.33333333));

    // Work around a bug in surfaceengine: All meshes must have cubes
    // or that engine can crash.
    Cube *cube = m_molecule->addCube();
    cube->setLimits(Eigen::Vector3d(0,0,0), Eigen::Vector3i(1,1,1), 1);
    std::vector<double> junkDVec (1);
    cube->setData(junkDVec);
    mesh->setCube(cube->id());

    return 0;
  }

  void TensorVis::setMolecule(Molecule *molecule)
  {
    m_molecule = molecule;
  }

  QDockWidget * TensorVis::dockWidget()
  {
    return 0;
  }

}

#include "tensorvis.moc"

Q_EXPORT_PLUGIN2(tensorvis, TensorVis::TensorVisFactory)