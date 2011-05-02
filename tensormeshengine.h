/**********************************************************************
  TensorMeshEngine - Engine for display meshes build from tensors

  Copyright (C) 2011 by David C. Lonie

  This source code is released under the New BSD License, (the "License").

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

 **********************************************************************/

#ifndef TENSORMESHENGINE_H
#define TENSORMESHENGINE_H

#include <avogadro/global.h>
#include <avogadro/engine.h>

namespace TensorVis {

  class TensorMeshEngine : public Avogadro::Engine
  {
      Q_OBJECT
      AVOGADRO_ENGINE("Tensor", tr("Tensor"),
                      tr("Renders tensor meshes built by the TensorVis extension."))

    public:
      TensorMeshEngine(QObject *parent=0);
      ~TensorMeshEngine();

      Engine *clone() const;
      bool renderOpaque(Avogadro::PainterDevice *pd);
      PrimitiveTypes primitiveTypes() const;

    private:
      bool renderPolygon(Avogadro::PainterDevice *pd,
                         Avogadro::Atom *a);
    };

  // Work around avogadro bug:
  using Avogadro::Plugin;

  class TensorMeshEngineFactory : public QObject,
    public Avogadro::PluginFactory
  {
    Q_OBJECT
    Q_INTERFACES(Avogadro::PluginFactory)
    AVOGADRO_ENGINE_FACTORY(TensorMeshEngine)
  };

}

#endif
