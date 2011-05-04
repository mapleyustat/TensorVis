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

#include "tensorviswidget.h"

#include <stdio.h>

namespace TensorVis
{

  TensorVisWidget::TensorVisWidget()
    : QDockWidget()
  {
    ui.setupUi(this);

    connect(ui.push_go, SIGNAL(clicked()),
            this, SIGNAL(generateMesh()));
    connect(this, SIGNAL(invalidInput()),
            this, SLOT(markAsInvalid()));
    connect(this, SIGNAL(validInput()),
            this, SLOT(markAsValid()));

    ui.edit_tensor->setCurrentFont(QFont("Monospace"));
    m_charFormat = ui.edit_tensor->textCursor().charFormat();

    m_origin << 0, 0, 0;
    m_tensor <<
      88.921500, 46.045900, -0.777100,
      17.744000, 18.222400, -0.248500,
      -0.900100, -0.432300, -81.653700;

    refreshGui();
  }

  TensorVisWidget::~TensorVisWidget()
  {
  }

  Eigen::Matrix3f TensorVisWidget::tensor()
  {
    return m_tensor;
  }

  Eigen::Vector3f TensorVisWidget::origin()
  {
    return m_origin;
  }

  Avogadro::Color3f TensorVisWidget::posColor()
  {
    float r=0.2, g=0.2, b=0.7;
    return Avogadro::Color3f(r,g,b);
  }

  Avogadro::Color3f TensorVisWidget::negColor()
  {
    float r=0.7, g=0.2, b=0.2;
    return Avogadro::Color3f(r,g,b);
  }

  void TensorVisWidget::refreshGui()
  {
    char text[256];
    snprintf(text, 256,
             "%9.5f %9.5f %9.5f\n"
             "%9.5f %9.5f %9.5f\n"
             "%9.5f %9.5f %9.5f\n",
             m_tensor(0, 0), m_tensor(0, 1), m_tensor(0, 2),
             m_tensor(1, 0), m_tensor(1, 1), m_tensor(1, 2),
             m_tensor(2, 0), m_tensor(2, 1), m_tensor(2, 2));

    ui.edit_tensor->blockSignals(true);
    ui.edit_tensor->setText(text);
    ui.edit_tensor->blockSignals(false);

    ui.edit_tensor->setCurrentFont(QFont("Monospace"));
  }

  void TensorVisWidget::validateEditor()
  {
    // Clear selection, otherwise there is a crash on Qt 4.7.2.
    QTextCursor tc = ui.edit_tensor->textCursor();
    tc.clearSelection();
    ui.edit_tensor->setTextCursor(tc);

    QString text = ui.edit_tensor->document()->toPlainText();
    QStringList lines = text.split("\n", QString::SkipEmptyParts);
    if (lines.size() != 3) {
      emit invalidInput();
      return;
    }

    const QRegExp parseRegExp
      ("\\s+|,|;|\\||\\[|\\]|\\{|\\}|\\(|\\)|\\&|/|<|>");

    QList<QStringList> stringVecs;
      for (int row = 0; row < 3; ++row) {
        stringVecs.append(lines.at(row).simplified()
                          .split(parseRegExp, QString::SkipEmptyParts));
        QStringList &stringVec = stringVecs[row];
        if (stringVec.size() != 3) {
          emit invalidInput();
          return;
        }
        for (int col = 0; col < 3; ++col) {
          bool ok;
          double val = stringVec[col].toDouble(&ok);
          if (!ok) {
            emit invalidInput();
            return;
          }
        }
      }

      emit validInput();
    }

  void TensorVisWidget::markAsInvalid()
  {
    QTextCursor tc (ui.edit_tensor->document());
    QTextCharFormat redFormat;
    redFormat.setForeground(QBrush(Qt::red));
    tc.movePosition(QTextCursor::Start);
    tc.movePosition(QTextCursor::End,
                    QTextCursor::KeepAnchor);
    ui.edit_tensor->blockSignals(true);
    tc.mergeCharFormat(redFormat);
    ui.edit_tensor->blockSignals(false);
    ui.edit_tensor->setCurrentFont(QFont("Monospace"));
  }

  void TensorVisWidget::markAsValid()
  {
    QTextCursor tc (ui.edit_tensor->document());
    tc.movePosition(QTextCursor::Start);
    tc.movePosition(QTextCursor::End,
                    QTextCursor::KeepAnchor);
    ui.edit_tensor->blockSignals(true);
    tc.setCharFormat(m_charFormat);
    ui.edit_tensor->blockSignals(false);
    ui.edit_tensor->setCurrentFont(QFont("Monospace"));
  }
}

#include "tensorviswidget.moc"
