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

#ifndef TENSORVIS_H
#define TENSORVIS_H

#include <avogadro/global.h>
#include <avogadro/plugin.h>
#include <avogadro/extension.h>

#include <Eigen/Core>

// This is a work around for a bug on older versions Avogadro, bug
// 3104853. Patch submitted.
using Avogadro::Plugin;

namespace TensorVis {

  class TensorVis : public Avogadro::Extension
  {
    Q_OBJECT
    AVOGADRO_EXTENSION("TensorVis",
                       tr("Plot a representation of a tensor over a molecule."),
                       tr("Plot a representation of a tensor over a molecule."))

  public:
    TensorVis(QObject *parent=0);
    ~TensorVis();

    virtual QList<QAction *> actions() const;
    virtual QString menuPath(QAction *action) const;
    virtual QUndoCommand* performAction(QAction *action,
                                        Avogadro::GLWidget *widget);
    virtual void setMolecule(Avogadro::Molecule *molecule);
    virtual QDockWidget * dockWidget();

  private:
    // "point" is on a sphere formed by iterating over phi and theta
    // in spherical coordinates:
    // x = r * sin(theta)*cos(phi)
    // y = r * sin(theta)*sin(phi)
    // z = r * cos(theta)
    static inline float calculateTensorValue(const Eigen::Vector3f &point,
                                             const Eigen::Matrix3f &tensor);

    QList<QAction *> m_actions;
    Avogadro::Molecule *m_molecule;
    QDockWidget *m_dock;

  };

  // Plugin factory setup
  class TensorVisFactory
    : public QObject,
      public Avogadro::PluginFactory
  {
    Q_OBJECT
    Q_INTERFACES(Avogadro::PluginFactory)
    AVOGADRO_EXTENSION_FACTORY(TensorVis)
  };

}

#endif

