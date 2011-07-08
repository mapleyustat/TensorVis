/**********************************************************************
  tensorviswidget.h Base class for crystal builder dockwidgets

  Copyright (C) 2011 by David C. Lonie

  This file is part of the Avogadro molecular editor project.
  For more information, see <http://avogadro.openmolecules.net/>

  This source code is released under the New BSD License, (the "License").

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
 ***********************************************************************/

#ifndef TENSORVISWIDGET_H
#define TENSORVISWIDGET_H

#include "../tensorvis.h"

#include <avogadro/color3f.h>

#include <Eigen/Core>

#include <QtGui/QTextCharFormat>
#include <QtGui/QDockWidget>

#include "ui_tensorviswidget.h"

namespace TensorVis
{
  class TensorVis;

  class TensorVisWidget : public QDockWidget
  {
    Q_OBJECT

  public:
    TensorVisWidget();
    virtual ~TensorVisWidget();

    Eigen::Matrix3f tensor();
    Eigen::Vector3f origin();
    Avogadro::Color3f posColor();
    Avogadro::Color3f negColor();

    const TensorVis::Resolution resolution() {
      return static_cast<TensorVis::Resolution>
          (this->ui.combo_resolution->currentIndex());}

    const double scale() {return this->ui.spin_scale->value();}

  signals:
    void clearMesh();
    void generateMesh();

  public slots:
    void refreshGui();

  protected slots:
    void useSelectedAtomAsOrigin();

  protected:
    Ui::TensorVisWidget ui;

  };

}

#endif
